//
// Created by 梦想家xixi on 2022/3/1.
//

#include <iomanip>
#include "VectorLoglessPairHMM.h"
#include "intel/pairhmm/IntelPairHmm.h"
#include "ReadUtils.h"
#include "haplotypecaller/ReadForPairHMM.h"
#include "parallel_hashmap/phmap.h"
#include "utils/pairhmm/PairHMMConcurrentControl.h"

VectorLoglessPairHMM::VectorLoglessPairHMM(PairHMMNativeArgumentCollection &args) : mHaplotypeDataArrayLength(0) {
	initNative(args.useDoublePrecision, args.pairHmmNativeThreads);
}

void VectorLoglessPairHMM::initialize(const std::vector<std::shared_ptr<Haplotype>> &haplotypes,
                                      const std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> &perSampleReadList,
                                      const int readMaxLength, const int haplotypeMaxLength) {

	mHaplotypeDataArrayLength = haplotypes.size();
	//mHaplotypeDataArray = std::shared_ptr<std::shared_ptr<HaplotypeDataHolder>[]>(new std::shared_ptr<HaplotypeDataHolder>[mHaplotypeDataArrayLength]);
	//mHaplotypeDataArray = std::shared_ptr<HaplotypeDataHolder[]>(new HaplotypeDataHolder[]());
	mHaplotypeDataArray.reserve(mHaplotypeDataArrayLength);
	int idx = 0;
	haplotypeToHaplotypeListIdxMap.clear();
	for (const std::shared_ptr<Haplotype> &currHaplotype: haplotypes) {
		mHaplotypeDataArray[idx] = HaplotypeDataHolder(currHaplotype->getBases().get(),
		                                               currHaplotype->getBasesLength());
		haplotypeToHaplotypeListIdxMap.emplace(currHaplotype, idx);
		idx++;
	}
}

void VectorLoglessPairHMM::computeLog10Likelihoods(SampleMatrix<SAMRecord, Haplotype> *logLikelihoods,
                                                   vector<shared_ptr<SAMRecord>> &processedReads,
                                                   phmap::flat_hash_map<SAMRecord *, shared_ptr<char[]>> *gcp) {
	if (processedReads.empty())
		return;

	int numReads = processedReads.size();

	// Get haplotypes
	std::vector<uint8_t *> haplotypes;
	std::vector<int> haplotypeLengths;
	haplotypes.reserve(mHaplotypeDataArrayLength);
	haplotypeLengths.reserve(mHaplotypeDataArrayLength);
	for (int i = 0; i < mHaplotypeDataArrayLength; ++i) {
		haplotypes.emplace_back(mHaplotypeDataArray[i].haplotypeBases);
		haplotypeLengths.emplace_back(mHaplotypeDataArray[i].length);
	}

	// An array that stores unique testcases and where all testcases are mapped to this array
	// Therefore, the number of testcases to be computed is reduced, but the length of mloglikelihoodarray is unchanged
	std::vector<testcase> uniqueTestcases;
	/* For Debugging*/
//	std::vector<testcase> old_uniqueTestcases;
	std::vector<int> mapAlltoUnique;
	uniqueTestcases.reserve(numReads * mHaplotypeDataArrayLength);
	mapAlltoUnique.reserve(numReads * mHaplotypeDataArrayLength);

	// Where do all the testcases related to a read start
	phmap::flat_hash_map<std::shared_ptr<ReadForPairHMM>, int, ReadForPairHMMHash, ReadForPairHMMEqual> uniqueReadForPairHMM;
	uniqueReadForPairHMM.reserve(numReads);

	// Generate unique testcases
	int _rslen;
	char *gapConts;
	uint8_t *reads, *readQuals;
	std::vector<shared_ptr<uint8_t[]>> insGops, delGops;
	insGops.reserve(numReads);
	delGops.reserve(numReads);
	for (int r = 0; r < numReads; ++r) {
		_rslen = processedReads[r]->getLength();
		insGops.emplace_back(ReadUtils::getBaseInsertionQualities(processedReads[r], _rslen));
		delGops.emplace_back(ReadUtils::getBaseDeletionQualities(processedReads[r], _rslen));
		gapConts = (*gcp)[processedReads[r].get()].get();
		readQuals = processedReads[r]->getBaseQualitiesNoCopy().get();
		reads = processedReads[r]->getBasesNoCopy().get();

		std::shared_ptr<ReadForPairHMM> readOfTestcase = std::make_shared<ReadForPairHMM>(_rslen, readQuals,
		                                                                                  insGops[r].get(),
		                                                                                  delGops[r].get(), gapConts,
		                                                                                  reads);
		auto readIt = uniqueReadForPairHMM.find(readOfTestcase);
		if (BOOST_LIKELY(readIt == uniqueReadForPairHMM.end())) {
			// Push testcases into uniqueTestcases and mark the index where the first testcase appears
			uniqueReadForPairHMM.emplace(readOfTestcase, uniqueTestcases.size());
			for (int h = 0; h < mHaplotypeDataArrayLength; ++h) {
				mapAlltoUnique.emplace_back(uniqueTestcases.size());
				uniqueTestcases.emplace_back(haplotypeLengths[h], haplotypes[h], readOfTestcase);
				/* For Debugging*/
//				old_uniqueTestcases.emplace_back(haplotypeLengths[h], haplotypes[h], readOfTestcase);
			}
		} else {
			// No element needs to be pushed into uniqueTestcases
			int mapStart = readIt->second;
			for (int h = 0; h < mHaplotypeDataArrayLength; ++h) {
				mapAlltoUnique.emplace_back(mapStart++);
				/* For Debugging*/
//				old_uniqueTestcases.emplace_back(haplotypeLengths[h], haplotypes[h], readOfTestcase);
			}
		}
	}

//	PairHMMConcurrentControl::unique_reads += uniqueReadForPairHMM.size();
//	PairHMMConcurrentControl::all_reads += numReads;
//	PairHMMConcurrentControl::unique_cases += uniqueTestcases.size();
//	PairHMMConcurrentControl::all_cases += numReads * mHaplotypeDataArrayLength;
//	if (uniqueTestcases.size() != numReads * mHaplotypeDataArrayLength) {
//		std::cout << "==========================\n";
//		std::cout << "read:\t" << uniqueReadForPairHMM.size() << " / " << numReads << std::endl;
//		std::cout << "case:\t" << uniqueTestcases.size() << " / " << numReads * mHaplotypeDataArrayLength << std::endl;
//	}

	//---for debugging
/*    for(int i=0; i<testcases.size(); i++)
    {
        for(int r=0; r<numReads; r++)
        {
            for(int i=0; i<testcases[r].rslen; i++)
            {
                cerr << testcases[r].rs[i];
            }
            cerr << endl;
            for(int i=0; i<testcases[r].rslen; i++)
            {
                cerr << (int)testcases[r].q[i] << " ";
            }
            cerr << endl;
            for(int i=0; i<testcases[r].rslen; i++)
            {
                cerr << (int)testcases[r].i[i] << " ";
            }
            cerr << endl;
            for(int i=0; i<testcases[r].rslen; i++)
            {
                cerr << (int)testcases[r].d[i] << " ";
            }
            cerr << endl;
            for(int i=0; i<testcases[r].rslen; i++)
            {
                cerr << (int)testcases[r].c[i] << " ";
            }
            cerr << endl;
            for(int i=0; i<testcases[r].haplen; i++)
            {
                cerr << testcases[r].hap[i];
            }
            cerr << endl;
            cerr << endl;
        }
    }*/

	// initialize probaility arrays and distm
	for (const auto &[readForPairHMM, _]: uniqueReadForPairHMM)
		readForPairHMM->initializeFloatVector();

	// Compute
	std::vector<double> uniqueLogLikelihoodArray(uniqueTestcases.size());
	computeLikelihoodsNative_concurrent(uniqueTestcases, uniqueLogLikelihoodArray);

	// Mapping results
	mLogLikelihoodArray.resize(mapAlltoUnique.size());
	for (int i = 0; i < mapAlltoUnique.size(); ++i)
		mLogLikelihoodArray[i] = uniqueLogLikelihoodArray[mapAlltoUnique[i]];

	/* For Debugging*/
//	std::vector<double> old_mLogLikelihoodArray(old_uniqueTestcases.size());
//	computeLikelihoodsNative_concurrent(old_uniqueTestcases, old_mLogLikelihoodArray);
//	assert(old_mLogLikelihoodArray == mLogLikelihoodArray);

	//---print the likelihoods calculated by PairHMM algorithm
//	cout << numReads << " " << mHaplotypeDataArrayLength << endl;
//	for (int r = 0; r < numReads; ++r) {
//		for (int h = 0; h < mHaplotypeDataArrayLength; ++h) {
//			cout.setf(ios::fixed);
//			cout << setprecision(5) << mLogLikelihoodArray[r * mHaplotypeDataArrayLength + h] << " ";
//		}
//		cout << endl;
//	}

	int readIdx = 0;
	for (int r = 0; r < numReads; ++r) {
		int hapIdx = 0;
		for (auto &haplotype: logLikelihoods->alleles()) {
			//Since the order of haplotypes in the List<Haplotype> and alleleHaplotypeMap is different,
			//get idx of current haplotype in the list and use this idx to get the right likelihoodValue
			int idxInsideHaplotypeList = haplotypeToHaplotypeListIdxMap.at(haplotype);
			logLikelihoods->set(hapIdx, r, mLogLikelihoodArray[readIdx + idxInsideHaplotypeList]);
			hapIdx++;
		}
		readIdx += mHaplotypeDataArrayLength;
	}
}