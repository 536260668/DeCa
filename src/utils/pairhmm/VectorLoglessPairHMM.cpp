//
// Created by 梦想家xixi on 2022/3/1.
//

#include "VectorLoglessPairHMM.h"
#include "intel/pairhmm/IntelPairHmm.h"
#include "ReadUtils.h"

VectorLoglessPairHMM::VectorLoglessPairHMM(PairHMMNativeArgumentCollection& args)
{
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
    for(const std::shared_ptr<Haplotype> & currHaplotype : haplotypes) {
        mHaplotypeDataArray[idx] = HaplotypeDataHolder(currHaplotype->getBases().get(), currHaplotype->getBasesLength());
        haplotypeToHaplotypeListIdxMap.insert(pair<std::shared_ptr<Haplotype>, int>(currHaplotype, idx));
        idx++;
    }
}

void VectorLoglessPairHMM::computeLog10Likelihoods(SampleMatrix<SAMRecord, Haplotype>* logLikelihoods,
                                              vector<shared_ptr<SAMRecord>>& processedReads,
                                              unordered_map<SAMRecord*, shared_ptr<char[]>>* gcp){
    if(processedReads.empty())
        return;

    int numReads = processedReads.size();


    std::vector<testcase> testcases;
    std::vector<uint8_t *> haplotypes;
    std::vector<int> haplotypeLengths;
    vector<shared_ptr<uint8_t[]>> insGops;
    vector<shared_ptr<uint8_t[]>> delGops;
    testcases.reserve(numReads * mHaplotypeDataArrayLength);

    // get haplotypes
    for(int i=0; i<mHaplotypeDataArrayLength; i++)
    {
        haplotypes.push_back(mHaplotypeDataArray[i].haplotypeBases);
        haplotypeLengths.push_back(mHaplotypeDataArray[i].length);
    }

    // get reads and create testcases
    for(int r=0; r<numReads; r++)
    {
        int length = processedReads[r]->getLength();
        uint8_t *reads = processedReads[r]->getBasesNoCopy().get();
        insGops.emplace_back((ReadUtils::getBaseInsertionQualities(processedReads[r], length)));
        delGops.emplace_back(ReadUtils::getBaseDeletionQualities(processedReads[r], length));
        char * gapConts = (*gcp)[processedReads[r].get()].get();
        uint8_t * readQuals = processedReads[r]->getBaseQualitiesNoCopy().get();

        for (int h = 0; h < mHaplotypeDataArrayLength; h++) {
            testcases.emplace_back(testcase(length, haplotypeLengths[h], readQuals, insGops[r].get(), delGops[r].get(), gapConts, haplotypes[h], reads));
        }
    }

    //---for debugging
/*
    for(int r=0; r<numReads; r++)
    {
        for(int i=0; i<testcases[r].rslen; i++)
        {
            cout << testcases[r].rs[i];
        }
        cout << endl;
        for(int i=0; i<testcases[r].rslen; i++)
        {
            cout << (int)testcases[r].q[i];
        }
        cout << endl;
        for(int i=0; i<testcases[r].rslen; i++)
        {
            cout << (int)testcases[r].i[i];
        }
        cout << endl;
        for(int i=0; i<testcases[r].rslen; i++)
        {
            cout << (int)testcases[r].d[i];
        }
        cout << endl;
        for(int i=0; i<testcases[r].rslen; i++)
        {
            cout << (int)testcases[r].c[i];
        }
        cout << endl;
        for(int i=0; i<testcases[r].haplen; i++)
        {
            cout << testcases[r].hap[i];
        }
        cout << endl;
    }
*/


    // TODO: finish this method 2022.5.3
    mLogLikelihoodArray.clear();
    mLogLikelihoodArray.reserve(numReads * mHaplotypeDataArrayLength);
    computeLikelihoodsNative(testcases, mLogLikelihoodArray);


}