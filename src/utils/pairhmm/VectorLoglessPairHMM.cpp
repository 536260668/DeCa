//
// Created by 梦想家xixi on 2022/3/1.
//

#include <iomanip>
#include "VectorLoglessPairHMM.h"
#include "intel/pairhmm/IntelPairHmm.h"
#include "ReadUtils.h"

VectorLoglessPairHMM::VectorLoglessPairHMM(PairHMMNativeArgumentCollection& args): mHaplotypeDataArrayLength(0)
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

    mLogLikelihoodArray.clear();
    mLogLikelihoodArray.reserve(numReads * mHaplotypeDataArrayLength);
    computeLikelihoodsNative(testcases, mLogLikelihoodArray);

    //---print the likelihoods calculated by PairHMM algorithm
/*    for(double & likelihood: mLogLikelihoodArray)
    {
        cerr.setf(ios::fixed);
        cerr << setprecision(5) << likelihood << " ";
    }
    cerr << endl;*/

    int readIdx = 0;
    for(int r=0; r<numReads; r++)
    {
        int hapIdx = 0;
        for(auto & haplotype: logLikelihoods->alleles())
        {
            //Since the order of haplotypes in the List<Haplotype> and alleleHaplotypeMap is different,
            //get idx of current haplotype in the list and use this idx to get the right likelihoodValue
            int idxInsideHaplotypeList = haplotypeToHaplotypeListIdxMap.at(haplotype);
            logLikelihoods->set(hapIdx, r, mLogLikelihoodArray[readIdx + idxInsideHaplotypeList]);
            hapIdx++;
        }
        readIdx += mHaplotypeDataArrayLength;
    }

}