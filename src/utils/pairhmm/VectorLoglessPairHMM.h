//
// Created by 梦想家xixi on 2022/3/1.
//

#ifndef MUTECT2CPP_MASTER_VECTORLOGLESSPAIRHMM_H
#define MUTECT2CPP_MASTER_VECTORLOGLESSPAIRHMM_H

#include "PairHMM.h"
#include "haplotype/Haplotype.h"
#include "samtools/SAMRecord.h"
#include "HaplotypeDataHolder.h"
#include "AssemblyResultSet.h"
#include "haplotypecaller/PairHMMNativeArgumentCollection.h"
#include "tiretree/buildTreeUtils.h"

class VectorLoglessPairHMM : public PairHMM{
private:
    //---two-dimensional array ?
    //std::shared_ptr<std::shared_ptr<HaplotypeDataHolder>[]> mHaplotypeDataArray;
    vector<HaplotypeDataHolder> mHaplotypeDataArray;
    phmap::flat_hash_map<std::shared_ptr<Haplotype>, int, hash_Haplotype, equal_Haplotype> haplotypeToHaplotypeListIdxMap;
    unsigned mHaplotypeDataArrayLength;
    vector<std::shared_ptr<Haplotype>> haps;
    tireTreeNode *root;

public:
    VectorLoglessPairHMM(PairHMMNativeArgumentCollection& args);

    /**
     * Create a VectorLoglessPairHMM
     *
     * @param implementation    which implementation to use (AVX or OMP)
     * @param args              arguments to the native GKL implementation
     */
    void initialize(const std::vector<std::shared_ptr<Haplotype>> & haplotypes,
                    const std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> & perSampleReadList,
                    int readMaxLength, int haplotypeMaxLength);

    void computeLog10Likelihoods(SampleMatrix<SAMRecord, Haplotype>* logLikelihoods,
                                 vector<shared_ptr<SAMRecord>>& processedReads,
                                 phmap::flat_hash_map<SAMRecord*, shared_ptr<char[]>>* gcp);

    void computeLog10Likelihoods_tiretree(SampleMatrix<SAMRecord, Haplotype>* logLikelihoods,
                                          vector<shared_ptr<SAMRecord>>& processedReads,
                                          phmap::flat_hash_map<SAMRecord*, shared_ptr<char[]>>* gcp);
};


#endif //MUTECT2CPP_MASTER_VECTORLOGLESSPAIRHMM_H
