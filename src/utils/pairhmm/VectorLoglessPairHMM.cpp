//
// Created by 梦想家xixi on 2022/3/1.
//

#include "VectorLoglessPairHMM.h"

void VectorLoglessPairHMM::initialize(const std::vector<std::shared_ptr<Haplotype>> &haplotypes,
                                      const std::map<std::string, std::vector<SAMRecord>> &perSampleReadList,
                                      const int readMaxLength, const int haplotypeMaxLength) {
    mHaplotypeDataArrayLength = haplotypes.size();
    mHaplotypeDataArray = std::shared_ptr<std::shared_ptr<HaplotypeDataHolder>[]>(new std::shared_ptr<HaplotypeDataHolder>[mHaplotypeDataArrayLength]);
    int idx = 0;
    haplotypeToHaplotypeListIdxMap.clear();
    for(const std::shared_ptr<Haplotype> & currHaplotype : haplotypes) {
        mHaplotypeDataArray[idx] = std::make_shared<HaplotypeDataHolder>(currHaplotype->getBases(), currHaplotype->getBasesLength());
    }
}
