//
// This is the main entry to Mutect2Cpp-master
// Created by lhh on 10/18/21.
//

#include <iostream>
#include <cstring>
#include <vector>
#include <queue>
#include "getopt.h"
#include "unistd.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "Mutect2Engine.h"
#include "QualityUtils.h"
#include "ReadFilter.h"
#include "M2ArgumentCollection.h"
#include "ActivityProfile.h"
#include "BandPassActivityProfile.h"
#include "variantcontext/GenotypeLikelihoods.h"
#include "samtools/SAMTools_decode.h"
#include "samtools/SamFileHeaderMerger.h"
#include "engine/ReferenceContext.h"
#include "read/ReadCache.h"
#include "MathUtils.h"
#include "intel/smithwaterman/IntelSmithWaterman.h"


static int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: Mutect2 -R reference.fa -I tumor.bam -I normal.bam -O output.vcf\n");
    fprintf(stderr, "Required Arguments:\n");
    fprintf(stderr, "--input,-I:String             BAM/SAM/CRAM file containing reads  This argument must be specified at least once.\n");
    fprintf(stderr, "                              Required\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--output,-O:File              The output recalibration table file to create  Required.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--reference,-R:String         Reference sequence file  Required. \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Optional Arguments:\n");
    fprintf(stderr, "--callable-depth              Minimum depth to be considered callable for Mutect stats. Does not affect genotyping. default: 10\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--max-prob-propagation-distance    Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--active-probability-threshold     Minimum probability for a locus to be considered active. default: 0.002\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--assembly-region-padding     Number of additional bases of context to include around each assembly region. default: 100\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--max-assembly-region-size    Maximum size of an assembly region. default: 300\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--min-assembly-region-size    Minimum size of an assembly region. default: 50\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-normal                        BAM sample name of normal(s), if any. May be URL-encoded as output by GetSampleName with -encode argument.\n");
    fprintf(stderr, "\n");
                    return EXIT_FAILURE;
}

// TODO: you can insert the filter in this method
// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1)
    {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
//        std::cout << " org : ";
//        for(int j = 0; j < b->core.l_qseq; j++) {
//            std::cout << static_cast<char>(bam_get_qual(b)[j] + '!');
//        }
//        std::cout << std::endl;
//        ReadFilter filter(b, aux->header);
        if ( ret<0 ) break;
//        if ( !filter.test() ) continue;
//        std::cout <<" flt : ";
//        for(int j = 0; j < b->core.l_qseq; j++) {
//            std::cout << static_cast<char>(bam_get_qual(b)[j] + '!');
//        }
//        std::cout << std::endl;
        break;
    }
    return ret;
}

int main(int argc, char *argv[])
{
    CigarOperatorUtils::initial();
    int c, n=0, reg_tid, tid;
    hts_pos_t beg, end, pos = -1, last_pos = -1;
    char * reg;
    char *in = nullptr, *output = nullptr, *ref = nullptr;
    std::vector<char*> input_bam;
    aux_t **data;
    bam_mplp_t mplp;
    uint32_t flags = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);

    int * n_plp;
    const bam_pileup1_t **plp;
    int ret;

    int baseQ = 1;
    int read_num = 0;

    M2ArgumentCollection MTAC = {10, 50, 0.002, 100, 50, 300};

    static struct option loptions[] =
            {
            {"input",       required_argument, NULL, 'I'},
            {"output",      required_argument, NULL, 'O'},
            {"reference",   required_argument, NULL, 'R'},
            {"callable-depth", required_argument, NULL, 1000},
            {"max-prob-propagation-distance", required_argument, NULL, 1001},
            {"active-probability-threshold", required_argument, NULL, 1002},
            {"assembly-region-padding", required_argument, NULL, 1003},
            {"max-assembly-region-size", required_argument, NULL, 1004},
            {"min-assembly-region-size", required_argument, NULL, 1005},
            {"normal", required_argument, NULL, 1006},
            { NULL, 0, NULL, 0 }
            };

    if (argc == 1 && isatty(STDIN_FILENO))
        return usage();

    while((c = getopt_long(argc, argv, "I:O:R:r:", loptions, NULL)) >= 0){
        switch (c) {
            case 'I':
                input_bam.emplace_back(strdup(optarg));
                n++;
                break;
            case 'O':
                output = strdup(optarg);
                break;
            case 'R':
                ref = strdup(optarg);
                break;
            case 'r':
                reg = strdup(optarg); break;   // parsing a region requires a BAM header
            case 1000:  //--callable-depth
                MTAC.callableDepth = atoi(optarg);
                break;
            case 1001:
                MTAC.maxProbPropagationDistance = atoi(optarg);
                break;
            case 1002:
                MTAC.activeProbThreshold = atof(optarg);
                break;
            case 1003:
                MTAC.assemblyRegionPadding = atoi(optarg);
                break;
            case 1004:
                MTAC.maxAssemblyRegionSize = atoi(optarg);
                break;
            case 1005:
                MTAC.minAssemblyRegionSize = atoi(optarg);
                break;
            case 1006:
                MTAC.normalSamples.insert(string(optarg));
                break;
        }
    }

    data = static_cast<aux_t **>(calloc(n, sizeof(aux_t *))); // data[i] for the i-th input
    reg_tid = 0; beg = 0; end = HTS_POS_MAX;  // set the default region
    vector<sam_hdr_t *> headers;// used to contain headers

    for (int i = 0; i < n; ++i) {
        data[i] = static_cast<aux_t *>(calloc(1, sizeof(aux_t)));
        data[i]->fp = hts_open(input_bam[i], "r"); // open BAM
        if (data[i]->fp == nullptr) {
            throw "Could not open sam/bam/cram files";
            exit(0);
        }
        data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        headers.emplace_back(data[i]->hdr);
        data[i]->flags = flags;
        data[i]->min_mapQ = Mutect2Engine::READ_QUALITY_FILTER_THRESHOLD;
    }
    SAMFileHeader* header = SAMTools_decode::merge_samFileHeaders(headers);
    for (int i = 0; i < n; ++i) {
        data[i]->header = header;
    }
    faidx_t * refPoint = fai_load3_format(ref, nullptr, nullptr, FAI_CREATE, FAI_FASTA);

    sam_hdr_t *h = data[0]->hdr; // easy access to the header of the 1st BAM   //---why use header of the first bam
    int nref = sam_hdr_nref(h);

    smithwaterman_initial();
    QualityUtils::initial();
    Mutect2Engine m2Engine(MTAC, ref, header);
    queue<std::shared_ptr<AssemblyRegion>> pendingRegions;
    ActivityProfile * activityProfile = new BandPassActivityProfile(MTAC.maxProbPropagationDistance, MTAC.activeProbThreshold, BandPassActivityProfile::MAX_FILTER_SIZE, BandPassActivityProfile::DEFAULT_SIGMA,
                                                                    true , header);
    int count = 0;
    // TODO: add multi-thread mode here
    for(int k=0; k<nref; k++)
    {
        std::string region = header->getSequenceDictionary().getSequences()[k].getSequenceName() + ":0-1999999";
        ReadCache cache(data, input_bam, k, region);
        int currentPose = 0;
        int len = 0;
        std::string contig = header->getSequenceDictionary().getSequences()[k].getSequenceName();
        char* refBases = faidx_fetch_seq(refPoint, contig.c_str(), 0, header->getSequenceDictionary().getSequences()[k].getSequenceLength(), &len);
        if(k != 0) {
            m2Engine.refCache.setTid(k);
        }
        while(cache.hasNextPos()) {
            AlignmentContext pileup = cache.getAlignmentContext();
            if(!activityProfile->isEmpty()){
                bool forceConversion = pileup.getPosition() != activityProfile->getEnd() + 1;
                vector<std::shared_ptr<AssemblyRegion>> * ReadyAssemblyRegions = activityProfile->popReadyAssemblyRegions(MTAC.assemblyRegionPadding, MTAC.minAssemblyRegionSize, MTAC.maxAssemblyRegionSize, forceConversion);
                for(const std::shared_ptr<AssemblyRegion>& newRegion : *ReadyAssemblyRegions)
                {
                    pendingRegions.emplace(newRegion);
                }
            }
//            if(pileup.getPosition() == 1000001)
//                std::cout << "hello";
            if(pileup.isEmpty()) {
                ActivityProfileState state(contig.c_str(), pileup.getPosition(), 0.0);
                activityProfile->add(state);
                continue;
            }
            SimpleInterval pileupInterval = SimpleInterval(contig, (int)pileup.getPosition(), (int)pileup.getPosition());
            ReferenceContext pileupRefContext(refBases, pileupInterval);

            ActivityProfileState profile = m2Engine.isActive(pileup, pileupRefContext);
            activityProfile->add(profile);

            if(!pendingRegions.empty() && IntervalUtils::isAfter(pileup.getLocation(), pendingRegions.front()->getExtendedSpan(), header->getSequenceDictionary())) {
                count++;
                if(count > 2000) {
                    break;
                }
                std::shared_ptr<AssemblyRegion> nextRegion = pendingRegions.front();
                pendingRegions.pop();
                Mutect2Engine::fillNextAssemblyRegionWithReads(nextRegion, cache);
                std::vector<std::shared_ptr<VariantContext>> variant = m2Engine.callRegion(nextRegion, pileupRefContext);
                if(variant.size() != 0)
                    std::cout << variant.size();
            }
        }


            // gather AlignmentContext to AssemblyRegion


        break;
    }

    //std::cout << "read_num : " << read_num << std::endl;


    // free the space


    return 0;
}