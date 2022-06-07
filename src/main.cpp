//
// This is the main entry to Mutect2Cpp-master
// Created by lhh on 10/18/21.
//

#include <iostream>
#include <cstring>
#include <thread>
#include <atomic>
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
#include "utils/BaseUtils.h"



struct Region{
    int _start;
    int _end;
	int _k;

	Region() {
		_start = 0;
		_end = 0;
		_k = 0;
	}

	Region(int start, int end, int k)
    {
        _start = start;
        _end = end;
		_k = k;
    }

	int getK() const {
		return _k;
	}

    int getStart() const
    {
        return _start;
    }

    int getEnd() const
    {
        return _end;
    }
};

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
	fprintf(stderr, "-T:Int                         Size of thread pool\n");
    fprintf(stderr, "\n");
	fprintf(stderr, "-L:String                      Specifies the name of the chromosome to be processed\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-M:String                      ML model path\n");
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

struct Shared{
	std::vector<std::shared_ptr<ReferenceCache>> refCaches;
	std::vector<char*> input_bam;
	SAMFileHeader* header{};
	std::string chromosomeName;
	M2ArgumentCollection MTAC;
	std::string modelPath;
	bool bqsr_within_mutect = false;
	std::shared_ptr<BQSRReadTransformer> tumorTransformer = nullptr;
	std::shared_ptr<BQSRReadTransformer> normalTransformer = nullptr;

	bool startFlag = false;
	std::vector<Region> regions;
	std::atomic<int> indexOfRegion{0};
	std::atomic<int> activeRegioncount{0};
	std::atomic<int> indexOfRefCache{0};
	std::atomic<int> numOfRefCache{0};
};

void threadFunc(Shared *w, char *ref, int n, int nref) {
	std::queue<std::shared_ptr<AssemblyRegion>> pendingRegions;
	ActivityProfile *activityProfile = new BandPassActivityProfile(w->MTAC.maxProbPropagationDistance, w->MTAC.activeProbThreshold, BandPassActivityProfile::MAX_FILTER_SIZE, BandPassActivityProfile::DEFAULT_SIGMA,true , w->header);
	VaraintAnnotatiorEngine annotatiorEngine;   // TODO: make it more elegant
	Mutect2Engine m2Engine(w->MTAC, w->header, w->modelPath, annotatiorEngine);
	std::vector<SAMSequenceRecord> headerSequences = w->header->getSequenceDictionary().getSequences();

	std::cout << "thread start\n";

	// create all refCaches
	// This operation is thread safe
	for (int k = w->indexOfRefCache++; k < nref; k = w->indexOfRefCache++) {
		if (!w->chromosomeName.empty() && headerSequences[k].getSequenceName() != w->chromosomeName)
			continue;
		//std::cout << "creating refCache \t" + headerSequences[k].getSequenceName() + "\n";
		w->refCaches[k] = std::make_shared<ReferenceCache>(ref, w->header, k);
		w->numOfRefCache++;
	}

	// Regenerate data to avoid conflicts
	auto data = static_cast<aux_t **>(calloc(n, sizeof(aux_t *))); // data[i] for the i-th input
	for (int i = 0; i < n; ++i) {
		data[i] = static_cast<aux_t *>(calloc(1, sizeof(aux_t)));
		data[i]->fp = hts_open(w->input_bam[i], "r"); // open BAM
		if (data[i]->fp == nullptr) {
			throw std::runtime_error("Could not open sam/bam/cram files");
		}
		data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
		data[i]->header = w->header;
	}

	// make sure all refCaches have been created
	while (BOOST_UNLIKELY(!w->startFlag)) std::this_thread::yield();
	if (BOOST_UNLIKELY(!w->chromosomeName.empty())) {
		while (w->numOfRefCache != 1) std::this_thread::yield();
	}
	else {
		while (w->numOfRefCache != nref) std::this_thread::yield();
	}

	Region tmpRegion;
	int start, end, k;
	for (int currentTask = w->indexOfRegion++; currentTask < w->regions.size(); currentTask = w->indexOfRegion++) {
		// running the task
		tmpRegion = w->regions[currentTask];
		start = tmpRegion.getStart();
		end =  tmpRegion.getEnd();
		k = tmpRegion.getK();
		std::string contig = headerSequences[k].getSequenceName();
		std::cout << "solving " + std::to_string(currentTask) + ' ' + contig + ' ' + std::to_string(start) + ' ' + std::to_string(end) + '\n';

		ReadCache cache(data, w->input_bam, k, start, end, w->MTAC.maxAssemblyRegionSize, w->refCaches[k], w->bqsr_within_mutect, w->tumorTransformer.get(), w->normalTransformer.get());
		m2Engine.setReferenceCache(w->refCaches[k].get());

		while(cache.hasNextPos()) {
			AlignmentContext pileup = cache.getAlignmentContext();
			/*if (pileup.getReadNum() != 0 && pileup.getPosition() >= start && pileup.getPosition() <= end){
				std::cout<<pileup.getLocation().getStart()+1<<" "<<pileup.getReadNum()<<std::endl;
			}*/
			if(!activityProfile->isEmpty()){
				bool forceConversion = pileup.getPosition() != activityProfile->getEnd() + 1;
				vector<std::shared_ptr<AssemblyRegion>> * ReadyAssemblyRegions = activityProfile->popReadyAssemblyRegions(w->MTAC.assemblyRegionPadding, w->MTAC.minAssemblyRegionSize, w->MTAC.maxAssemblyRegionSize, forceConversion);
				for(const std::shared_ptr<AssemblyRegion>& newRegion : *ReadyAssemblyRegions)
				{
					if(newRegion->getStart() >= start && newRegion->getStart() < end)
						pendingRegions.emplace(newRegion);
				}
			}

			if(pileup.isEmpty()) {
				std::shared_ptr<ActivityProfileState> state = std::make_shared<ActivityProfileState>(contig.c_str(), pileup.getPosition(), 0.0);
				activityProfile->add(state);
				continue;
			}
			std::shared_ptr<SimpleInterval> pileupInterval = std::make_shared<SimpleInterval>(contig, (int)pileup.getPosition(), (int)pileup.getPosition());
			char refBase = w->refCaches[k]->getBase(pileup.getPosition());
			ReferenceContext pileupRefContext(pileupInterval, refBase); //---is this variable necessary?

			std::shared_ptr<ActivityProfileState> profile = m2Engine.isActive(pileup);
			activityProfile->add(profile);

			if(!pendingRegions.empty() && IntervalUtils::isAfter(pileup.getLocation(), *pendingRegions.front()->getExtendedSpan(), w->header->getSequenceDictionary())) {
				w->activeRegioncount++;

				std::shared_ptr<AssemblyRegion> nextRegion = pendingRegions.front();

				//---print the region
				//std::cout << nextRegion->getContig() + " " + to_string(nextRegion->getStart()+1) + " " + to_string(nextRegion->getEnd()+1) + '\n';
				pendingRegions.pop();

				Mutect2Engine::fillNextAssemblyRegionWithReads(nextRegion, cache);
				std::vector<std::shared_ptr<VariantContext>> variant = m2Engine.callRegion(nextRegion, pileupRefContext);
			}
		}

		// pop the AssemblyRegion left
		while(!activityProfile->isEmpty())
		{
			vector<std::shared_ptr<AssemblyRegion>> * ReadyAssemblyRegions = activityProfile->popReadyAssemblyRegions(w->MTAC.assemblyRegionPadding, w->MTAC.minAssemblyRegionSize, w->MTAC.maxAssemblyRegionSize, true);
			for(const std::shared_ptr<AssemblyRegion>& newRegion : *ReadyAssemblyRegions)
			{
				if(newRegion->getStart() >= start && newRegion->getStart() < end)
					pendingRegions.emplace(newRegion);
			}
		}

		while(!pendingRegions.empty())
		{
			w->activeRegioncount++;
			std::shared_ptr<AssemblyRegion> nextRegion = pendingRegions.front();

			//---print the region
			//std::cout << nextRegion->getContig() + " " + to_string(nextRegion->getStart()+1) + " " + to_string(nextRegion->getEnd()+1) + '\n';

			pendingRegions.pop();
			Mutect2Engine::fillNextAssemblyRegionWithReads(nextRegion, cache);
			// ReferenceContext is not needed for the time being
			std::shared_ptr<SimpleInterval> pileupInterval = std::make_shared<SimpleInterval>(contig, 0, 0);
			ReferenceContext tmp{pileupInterval, N};
			std::vector<std::shared_ptr<VariantContext>> variant = m2Engine.callRegion(nextRegion, tmp); // TODO: callRegion() needs pileupRefContext
		}
	}

	// free the space
	delete activityProfile;
	for(int i = 0; i < n; i++)
	{
		hts_close(data[i]->fp);
		sam_hdr_destroy(data[i]->hdr);
		free(data[i]);
	}
	free(data);
	std::cout << "thread exit\n";
}

int main(int argc, char *argv[])
{
    CigarOperatorUtils::initial();
    int c, n = 0;
	char * reg;
    char *output = nullptr, *ref = nullptr;
	aux_t **data;
	Shared sharedData;
    char * tumor_table = nullptr, * normal_table = nullptr;
	int thread_num = 1;

	sharedData.MTAC = {10, 50, 0.002, 100, 50, 300, ""};

    static struct option loptions[] =
            {
            {"input",       required_argument, nullptr, 'I'},
            {"output",      required_argument, nullptr, 'O'},
            {"reference",   required_argument, nullptr, 'R'},
            {"thread",      required_argument, nullptr, 'T'},
            {"chromosome",  required_argument, nullptr, 'L'},
            {"model",       required_argument, nullptr, 'M'},
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

    while((c = getopt_long(argc, argv, "I:O:R:r:T:L:M:", loptions, nullptr)) >= 0){
        switch (c) {
            case 'I':
	            sharedData.input_bam.emplace_back(strdup(optarg));
                n++;
                break;
            case 'O':
                output = strdup(optarg);
                break;
            case 'R':
                ref = strdup(optarg);
                break;
	        case 'T':
				thread_num = std::max(1, atoi(optarg));
		        break;
	        case 'L':
		        sharedData.chromosomeName = std::string (optarg);
		        break;
	        case 'M':
		        sharedData.modelPath = std::string (optarg);
		        break;
            case 'r':
                reg = strdup(optarg);
                break;   // parsing a region requires a BAM header
            case 1000:  //--callable-depth
	            sharedData.MTAC.callableDepth = atoi(optarg);
                break;
            case 1001:
	            sharedData.MTAC.maxProbPropagationDistance = atoi(optarg);
                break;
            case 1002:
	            sharedData.MTAC.activeProbThreshold = atof(optarg);
                break;
			case 1003:
	            sharedData.MTAC.assemblyRegionPadding = atoi(optarg);
                break;
            case 1004:
	            sharedData.MTAC.maxAssemblyRegionSize = atoi(optarg);
                break;
            case 1005:
	            sharedData.MTAC.minAssemblyRegionSize = atoi(optarg);
                break;
            case 1006:
	            sharedData.MTAC.normalSample = string(optarg);
                break;
            case 1007:
                sharedData.bqsr_within_mutect = true;
                break;
            case 1008:
                tumor_table = strdup(optarg);
                break;
            case 1009:
                normal_table = strdup(optarg);
                break;
        }
    }

    adjust_input_bam(sharedData.input_bam, sharedData.MTAC.normalSample);

	data = static_cast<aux_t **>(calloc(n, sizeof(aux_t *))); // data[i] for the i-th input
    vector<sam_hdr_t *> headers;// used to contain headers

    for (int i = 0; i < n; ++i) {
	    data[i] = static_cast<aux_t *>(calloc(1, sizeof(aux_t)));
	    data[i]->fp = hts_open(sharedData.input_bam[i], "r"); // open BAM
        if (data[i]->fp == nullptr) {
            throw std::runtime_error("Could not open sam/bam/cram files");
        }
	    data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        headers.emplace_back(data[i]->hdr);
    }
    sharedData.header = SAMTools_decode::merge_samFileHeaders(headers);
    for (int i = 0; i < n; ++i) {
	    data[i]->header = sharedData.header;
    }
    faidx_t * refPoint = fai_load3_format(ref, nullptr, nullptr, FAI_CREATE, FAI_FASTA);

    sam_hdr_t *h = data[0]->hdr; // easy access to the header of the 1st BAM   //---why use header of the first bam
    int nref = sam_hdr_nref(h);

    smithwaterman_initial();
    QualityUtils::initial();
	BaseUtils::initial();

    if(sharedData.bqsr_within_mutect)
    {
        ApplyBQSRArgumentCollection bqsrArgs;
	    sharedData.tumorTransformer = std::make_shared<BQSRReadTransformer>(tumor_table, bqsrArgs);
	    sharedData.normalTransformer = std::make_shared<BQSRReadTransformer>(normal_table, bqsrArgs);
    }

	//start threads
	std::vector<std::thread> threads;
	sharedData.refCaches.resize(nref);
	for (int i = 0; i < thread_num; ++i) {
		threads.emplace_back(&threadFunc, &sharedData, ref, n, nref);
	}

	std::vector<SAMSequenceRecord> headerSequences = sharedData.header->getSequenceDictionary().getSequences();
    for(int k = 0; k < nref; k++)
    {
	    if (!sharedData.chromosomeName.empty() && headerSequences[k].getSequenceName() != sharedData.chromosomeName)
		    continue;
        hts_pos_t ref_len = sam_hdr_tid2len(data[0]->hdr, k);   // the length of reference sequence
        int start = 0;
        int end = ref_len < REGION_SIZE - 1 ? ref_len : REGION_SIZE - 1;
        while(end != ref_len)
        {
	        sharedData.regions.emplace_back(start, end, k);
            start += REGION_SIZE;
            end = end + REGION_SIZE < ref_len ? end + REGION_SIZE : ref_len;
        }
	    sharedData.regions.emplace_back(start, end, k);
    }

	std::cout << "count of regions: " + std::to_string(sharedData.regions.size()) + "\n";
	sharedData.startFlag = true;

	for(auto &thread : threads) {
		thread.join();
	}

    // free the space
    fai_destroy(refPoint);
    for(char * input_file : sharedData.input_bam)
        free(input_file);
    free(output);
    free(ref);
    free(tumor_table);
    free(normal_table);
    delete sharedData.header;
    for(int i = 0; i < n; i++)
    {
        hts_close(data[i]->fp);
        sam_hdr_destroy(data[i]->hdr);
        free(data[i]);
    }
    free(data);

    return 0;
}