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
#include "Mutect2Engine.h"
#include "QualityUtils.h"
#include "ReadFilter.h"
#include "AlignmentContext.h"
#include "M2ArgumentCollection.h"
#include "ActivityProfile.h"
#include "BandPassActivityProfile.h"
#include "variantcontext/GenotypeLikelihoods.h"

typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    sam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    int min_mapQ; // mapQ filter;
    uint32_t flags;  // read filtering flags
} aux_t;



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
        if ( ret<0 ) break;
        if ( (int)b->core.qual < aux->min_mapQ || (int)b->core.qual == QualityUtils::MAPPING_QUALITY_UNAVALIABLE) continue;
        if ( b->core.flag & aux->flags ) continue;
        if ( !ReadFilter::ReadLengthTest(b) ) continue;
        break;
    }
    return ret;
}

int main(int argc, char *argv[])
{

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
    vector<sam_hdr_t *> headers;  // used to contain headers

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



    sam_hdr_t *h = data[0]->hdr; // easy access to the header of the 1st BAM   //---why use header of the first bam
    int nref = sam_hdr_nref(h);

    Mutect2Engine m2Engine(MTAC, ref);
    ActivityProfile * activityProfile = new BandPassActivityProfile(MTAC.maxProbPropagationDistance, MTAC.activeProbThreshold, BandPassActivityProfile::MAX_FILTER_SIZE, BandPassActivityProfile::DEFAULT_SIGMA,
                                                                    true , h);
    queue<AssemblyRegion> pendingRegions;

    // TODO: add multi-thread mode here
    for(int k=0; k<nref; k++)
    {

        // set the iterator
        for (int i = 0; i < n; ++i) {
            hts_idx_t *idx = nullptr;
            idx = sam_index_load(data[i]->fp, input_bam[i]);
            //data[i]->iter = sam_itr_queryi(idx, k, 0, sam_hdr_tid2len(h, k)); // set the iterator
            data[i]->iter = sam_itr_queryi(idx, k, 0, 300000);  //---only for test
            hts_idx_destroy(idx); // the index is not needed any more; free the memory
        }

        // the core multi-pileup loop
        mplp = bam_mplp_init(n, read_bam, (void**)data);

        if ( bam_mplp_init_overlaps(mplp) < 0) {
            throw "failed to init overlap detection";
        }

        n_plp = static_cast<int *>(calloc(n, sizeof(int))); // n_plp[i] is the number of covering reads from the i-th BAM
        plp = static_cast<const bam_pileup1_t **>(calloc(n, sizeof(bam_pileup1_t *))); // plp[i] points to the array of covering reads (internal in mplp)
        // if there is no covered position, ret = 0 and pos = -1
        while ((ret=bam_mplp64_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position

            if (pos < beg || pos >= end) continue; // out of range; skip
            if (tid >= sam_hdr_nref(h)) continue;     // diff number of @SQ lines per file?


/*            //fputs(sam_hdr_tid2name(h, tid), stdout);
            //fprintf(stdout, "\t%" PRIhts_pos, pos+1); // a customized printf() would be faster
            for (int i = 0; i < n; ++i) { // base level filters have to go here
                int j, m = 0;
                for (j = 0; j < n_plp[i]; ++j) {
                    const bam_pileup1_t *p = plp[i] + j; // DON'T modify plp[][] unless you really know
                    if ((p->is_del) || p->is_refskip) ++m; // having dels or refskips at tid:pos
                    else if (p->qpos < p->b->core.l_qseq &&
                    bam_get_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality

                    //printf("%s\n", bam_aux2Z(bam_aux_get(plp[i]->b, "RG")));
                }
                //fprintf(stdout, "\t%d", n_plp[i] - m); // this the depth to output
                read_num += n_plp[i];
            }*/



            AlignmentContext pileup(headers, tid, pos, n, n_plp, plp);
            if(!activityProfile->isEmpty()){
                bool forceConversion = pileup.getPosition() != activityProfile->getEnd() + 1;
                vector<AssemblyRegion> * ReadyAssemblyRegions = activityProfile->popReadyAssemblyRegions(MTAC.assemblyRegionPadding, MTAC.minAssemblyRegionSize, MTAC.maxAssemblyRegionSize, forceConversion);
                for(AssemblyRegion & region : *ReadyAssemblyRegions)
                {
                    pendingRegions.emplace(region);
                }
            }

            ActivityProfileState profile = m2Engine.isActive(&pileup);

            activityProfile->add(profile);



            // gather AlignmentContext to AssemblyRegion
        }


    }

    //std::cout << "read_num : " << read_num << std::endl;


    // free the space
    delete activityProfile;
    for (char * file_name: input_bam) {
        free(file_name);
    }
    for (int i = 0; i < n && data[i]; ++i) {
        sam_hdr_destroy(data[i]->hdr);
        if (data[i]->fp) sam_close(data[i]->fp);
        hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data);


    return 0;
}