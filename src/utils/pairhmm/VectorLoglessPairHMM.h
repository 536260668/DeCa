//
// Created by 梦想家xixi on 2022/3/1.
//

#ifndef MUTECT2CPP_MASTER_VECTORLOGLESSPAIRHMM_H
#define MUTECT2CPP_MASTER_VECTORLOGLESSPAIRHMM_H

#include "haplotype/Haplotype.h"
#include "samtools/SAMRecord.h"
#include "HaplotypeDataHolder.h"
#include "AssemblyResultSet.h"

class VectorLoglessPairHMM {
private:
    std::shared_ptr<std::shared_ptr<HaplotypeDataHolder>[]> mHaplotypeDataArray;
    std::unordered_map<std::shared_ptr<Haplotype>, int, hash_Haplotype, equal_Haplotype> haplotypeToHaplotypeListIdxMap;
    unsigned mHaplotypeDataArrayLength;

public:
    /**
     * Create a VectorLoglessPairHMM
     *
     * @param implementation    which implementation to use (AVX or OMP)
     * @param args              arguments to the native GKL implementation
     */
    void initialize(const std::vector<std::shared_ptr<Haplotype>> & haplotypes,
                    const std::map<std::string, std::vector<SAMRecord>> & perSampleReadList,
                    int readMaxLength, int haplotypeMaxLength);
};


#endif //MUTECT2CPP_MASTER_VECTORLOGLESSPAIRHMM_H
