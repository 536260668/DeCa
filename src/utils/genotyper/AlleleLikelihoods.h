//
// Created by lhh on 4/23/22.
//

#ifndef MUTECT2CPP_MASTER_ALLELELIKELIHOODS_H
#define MUTECT2CPP_MASTER_ALLELELIKELIHOODS_H

#include <vector>
#include <string>
#include <map>
#include <cassert>
#include "Allele.h"
#include "samtools/SAMRecord.h"

using namespace std;

template <typename E, typename A>
class SampleMatrix;

template <typename E, typename A>
class AlleleLikelihoods {
    friend class SampleMatrix<E,A>;

private:
    const static int MISSING_REF = -1;

    /**
    * Index of the reference allele if any, otherwise {@link #MISSING_REF}.
    */
    int referenceAlleleIndex = MISSING_REF;

    /**
     * Sample matrices lazily initialized (the elements not the array) by invoking {@link #sampleMatrix(int)}.
     */
    vector<SampleMatrix<E, A>> sampleMatrices;

    // Search for the reference allele, if not found the index is {@link MISSING_REF}.
    static int findReferenceAllele(std::vector<shared_ptr<A>>& alleles){
        int number = alleles.size();
        for(int i=0; i<number; i++)
        {
            if(alleles[i]->getIsReference())
                return i;
        }
        return MISSING_REF;
    }

    void setupIndexes(std::map<std::string, std::vector<std::shared_ptr<E>>>& evidenceBySample, int sampleCount, int alleleCount) {
        for(int s = 0; s < sampleCount; s++)
        {
            std::string& sample = samples[s];

            evidenceBySampleIndex.emplace_back(evidenceBySample[sample]);

            int sampleEvidenceCount = evidenceBySampleIndex[s].size();

            valuesBySampleIndex.emplace_back(vector<vector<double>>(alleleCount));
        }
    }

protected:
    bool isNaturalLog = false;

    /**
     * Evidence by sample index. Each sub array contains reference to the evidence of the ith sample.
     */
    std::vector<std::vector<std::shared_ptr<E>>> evidenceBySampleIndex;

    std::vector<std::vector<std::vector<double>>> valuesBySampleIndex;

    std::vector<std::string>& samples;

    vector<shared_ptr<A>>& alleles;



public:
    const static double LOG_10_INFORMATIVE_THRESHOLD;
    const static double NATURAL_LOG_INFORMATIVE_THRESHOLD;

    AlleleLikelihoods(vector<string>& samples, vector<shared_ptr<A>>& alleles, map<string, vector<std::shared_ptr<SAMRecord>>>& evidenceBySample) :
         samples(samples), alleles(alleles)
    {
        int sampleCount = samples.size();
        int alleleCount = alleles.size();

        evidenceBySampleIndex.reserve(sampleCount);
        valuesBySampleIndex.reserve(sampleCount);
        referenceAlleleIndex = findReferenceAllele(alleles);

        setupIndexes(evidenceBySample, sampleCount, alleleCount);
        sampleMatrices.reserve(sampleCount);
        for(int i=0; i<sampleCount; i++)
        {
            sampleMatrices.template emplace_back(i, this);
        }
    }

    /**
     * Returns the units of evidence that belong to a sample sorted by their index (within that sample).
     *
     * @param sampleIndex the requested sample.
     * @return never {@code null} but perhaps a zero-length array if there is no evidence in sample. No element in
     *   the array will be null.
     */
    vector<shared_ptr<E>>& sampleEvidence(int sampleIndex){
        return evidenceBySampleIndex[sampleIndex];
    }


    /**
     * Returns an evidence vs allele likelihood matrix corresponding to a sample.
     */
    SampleMatrix<E, A> * sampleMatrix(int sampleIndex)
    {
        assert(sampleIndex >= 0 && sampleIndex < samples.size());

        return &sampleMatrices[sampleIndex];
    }


};

template <typename E, typename A>
class SampleMatrix{
private:
    int sampleIndex;
    AlleleLikelihoods<E, A> * likelihood;

public:
    SampleMatrix()
    {
        sampleIndex = 0;
        likelihood = nullptr;
    }

    SampleMatrix(int sampleIndex, AlleleLikelihoods<E, A> * likelihoods)
    {
        this->sampleIndex = sampleIndex;
        this->likelihood = likelihoods;
    }

    //std::vector<E> evidence()

    void set(AlleleLikelihoods<E, A> &a,int alleleIndex, int evidenceIndex, double value){
        a.valuesBySampleIndex[sampleIndex][alleleIndex][evidenceIndex] = value;
    }

    double get(AlleleLikelihoods<E, A> &a, int alleleIndex, int evidenceIndex){
        return a.valuesBySampleIndex[sampleIndex][alleleIndex][evidenceIndex];
    }

    vector<shared_ptr<E>> & evidence(){
        return likelihood->sampleEvidence(sampleIndex);
    }

};


#endif //MUTECT2CPP_MASTER_ALLELELIKELIHOODS_H
