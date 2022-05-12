//
// This is the main entry to Mutect2Cpp-master
// Created by lhh on 10/18/21.
//

#include <iostream>
#include <cstring>
#include <vector>
#include <queue>
#include <cassert>
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
#include "intel/smithwaterman/IntelSmithWaterman.h"
#include "ReferenceCache.h"


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
    fprintf(stderr, "--bqsr-with-mutect             Add ApplyBQSR into this tool. if set, --tumor-table and --normal-table is needed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--tumor-table                  Recalibration table for tumor reads, generated in the BaseRecalibrator part\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--normal-table                 Recalibration table for normal reads, generated in the BaseRecalibrator part\n");
    fprintf(stderr, "\n");
                    return EXIT_FAILURE;
}

/**
 * adjust input_bam in main() to make input_bam[0] normal bam
 * TODO: Maybe this method can be more elegant
 */
void adjust_input_bam(std::vector<char*> & input_bam, std::string & normalSample)
{
    for(int i=0; i<input_bam.size(); i++)
    {
        if(strstr(input_bam[i], normalSample.c_str()))
        {
            if(i != 0){
                char * temp = input_bam[0];
                input_bam[0] = input_bam[i];
                input_bam[i] = temp;
            }
            break;
        }
    }

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

    int read_num = 0;
    bool bqsr_within_mutect = false;    //---Maybe these parameters can make a struct
    char * tumor_table = nullptr, * normal_table = nullptr;

    M2ArgumentCollection MTAC = {10, 50, 0.002, 100, 50, 300, ""};

    static struct option loptions[] =
            {
            {"input",       required_argument, nullptr, 'I'},
            {"output",      required_argument, nullptr, 'O'},
            {"reference",   required_argument, nullptr, 'R'},
            {"callable-depth", required_argument, nullptr, 1000},
            {"max-prob-propagation-distance", required_argument, nullptr, 1001},
            {"active-probability-threshold", required_argument, nullptr, 1002},
            {"assembly-region-padding", required_argument, nullptr, 1003},
            {"max-assembly-region-size", required_argument, nullptr, 1004},
            {"min-assembly-region-size", required_argument, nullptr, 1005},
            {"normal", required_argument, nullptr, 1006},
            {"bqsr-within-mutect", optional_argument, nullptr, 1007},
            {"tumor-table", required_argument, nullptr, 1008},
            {"normal-table", required_argument, nullptr, 1009},
            { nullptr, 0, nullptr, 0 }
            };

    if (argc == 1 && isatty(STDIN_FILENO))
        return usage();

    while((c = getopt_long(argc, argv, "I:O:R:r:", loptions, nullptr)) >= 0){
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
                reg = strdup(optarg);
                break;   // parsing a region requires a BAM header
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
                MTAC.normalSample = string(optarg);
                assert(MTAC.normalSample != "");
                break;
            case 1007:
                bqsr_within_mutect = true;
                break;
            case 1008:
                tumor_table = strdup(optarg);
                break;
            case 1009:
                normal_table = strdup(optarg);
                break;
        }
    }


    adjust_input_bam(input_bam, MTAC.normalSample);

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
    ActivityProfile * activityProfile = new BandPassActivityProfile(MTAC.maxProbPropagationDistance, MTAC.activeProbThreshold, BandPassActivityProfile::MAX_FILTER_SIZE, BandPassActivityProfile::DEFAULT_SIGMA,true , header);
    std::shared_ptr<BQSRReadTransformer> tumorTransformer(nullptr);
    std::shared_ptr<BQSRReadTransformer> normalTransformer(nullptr);
    if(bqsr_within_mutect)
    {
        ApplyBQSRArgumentCollection bqsrArgs;
        tumorTransformer = std::make_shared<BQSRReadTransformer>(tumor_table, bqsrArgs);
        normalTransformer = std::make_shared<BQSRReadTransformer>(normal_table, bqsrArgs);
    }


    int count = 0;
    // TODO: add multi-thread mode here
    for(int k=18; k<nref; k++)
    {
        hts_pos_t ref_len = sam_hdr_tid2len(data[0]->hdr, k);   // the length of reference sequence

        int len = ref_len < REGION_SIZE ? ref_len : REGION_SIZE;
        std::string region = std::string(sam_hdr_tid2name(data[0]->hdr, k)) + ":0-" + to_string(len);
        std::string contig = std::string(sam_hdr_tid2name(data[0]->hdr, k));
        std::shared_ptr<ReferenceCache>  refCache = std::make_shared<ReferenceCache>(ref, data[0]->header, k);
        ReadCache cache(data, input_bam, k, region, refCache, bqsr_within_mutect, tumorTransformer.get(), normalTransformer.get());

        m2Engine.setReferenceCache(refCache.get());
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

            if(pileup.isEmpty()) {
                std::shared_ptr<ActivityProfileState> state = std::make_shared<ActivityProfileState>(contig.c_str(), pileup.getPosition(), 0.0);
                activityProfile->add(state);
                continue;
            }
            std::shared_ptr<SimpleInterval> pileupInterval = std::make_shared<SimpleInterval>(contig, (int)pileup.getPosition(), (int)pileup.getPosition());
            char refBase = refCache->getBase(pileup.getPosition());
            ReferenceContext pileupRefContext(pileupInterval, refBase);

            std::shared_ptr<ActivityProfileState> profile = m2Engine.isActive(pileup);
            activityProfile->add(profile);

            if(!pendingRegions.empty() && IntervalUtils::isAfter(pileup.getLocation(), *pendingRegions.front()->getExtendedSpan(), header->getSequenceDictionary())) {
                count++;

                std::shared_ptr<AssemblyRegion> nextRegion = pendingRegions.front();

/*                if(count > 20) {
                    //std::cout << *nextRegion;
                    break;
                }*/
                pendingRegions.pop();

                Mutect2Engine::fillNextAssemblyRegionWithReads(nextRegion, cache);
                std::vector<std::shared_ptr<VariantContext>> variant = m2Engine.callRegion(nextRegion, pileupRefContext);
                //if(variant.size() != 0)
                  //  std::cout << variant.size();
            }
        }

        // pop the AssemblyRegion left
        while(!activityProfile->isEmpty())
        {
            vector<std::shared_ptr<AssemblyRegion>> * ReadyAssemblyRegions = activityProfile->popReadyAssemblyRegions(MTAC.assemblyRegionPadding, MTAC.minAssemblyRegionSize, MTAC.maxAssemblyRegionSize, true);
            for(const std::shared_ptr<AssemblyRegion>& newRegion : *ReadyAssemblyRegions)
            {
                pendingRegions.emplace(newRegion);
            }
        }

        while(!pendingRegions.empty())
        {
            count++;
            std::shared_ptr<AssemblyRegion> nextRegion = pendingRegions.front();
            //std::cout << *nextRegion;
            pendingRegions.pop();

            Mutect2Engine::fillNextAssemblyRegionWithReads(nextRegion, cache);
            //std::vector<std::shared_ptr<VariantContext>> variant = m2Engine.callRegion(nextRegion, pileupRefContext); // TODO: callRegion() needs pileupRefContext
        }
        break;
    }

    //std::cout << "read_num : " << read_num << std::endl;


    // free the space
    delete activityProfile;
    fai_destroy(refPoint);
    for(char * input_file : input_bam)
        free(input_file);
    free(output);
    free(ref);
    free(tumor_table);
    free(normal_table);
    delete header;
    for(int i=0; i<n; i++)
    {
        hts_close(data[i]->fp);
        sam_hdr_destroy(data[i]->hdr);
        free(data[i]);
    }
    free(data);

    return 0;
}