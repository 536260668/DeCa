//
// Created by lhh on 4/23/22.
//

#ifndef MUTECT2CPP_MASTER_ALLELELIKELIHOODS_H
#define MUTECT2CPP_MASTER_ALLELELIKELIHOODS_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <cassert>
#include <functional>
#include <unordered_set>
#include "Allele.h"
#include "samtools/SAMRecord.h"
#include "MathUtils.h"

using namespace std;

template <typename E, typename A>
class SampleMatrix;
template <typename E, typename A>
class BestAllele;

template <typename E, typename A>
class AlleleLikelihoods {
    friend class SampleMatrix<E,A>;
    friend class BestAllele<E, A>;

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

    double getInformativeThreshold() {
        return isNaturalLog ? NATURAL_LOG_INFORMATIVE_THRESHOLD : LOG_10_INFORMATIVE_THRESHOLD;
    }

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

            valuesBySampleIndex.emplace_back(vector<vector<double>>(alleleCount, vector<double>(sampleEvidenceCount, 0.0)));
        }
    }



    /**
     * Search the best allele for a unit of evidence.
     *
     * @param sampleIndex including sample index.
     * @param evidenceIndex  target evidence index.
     *
     * @param priorities An array of allele priorities (higher values have higher priority) to be used, if present, to break ties for
     *                   uninformative likelihoods, in which case the evidence is assigned to the allele with the higher score.
     * @return never {@code null}, but with {@link BestAllele#allele allele} == {@code null}
     * if non-could be found.
     */
    shared_ptr<BestAllele<E, A>> searchBestAllele(int sampleIndex, int evidenceIndex, bool canBeReference, double* priorities){
        int alleleCount = alleles->size();
        if (alleleCount == 0 || (alleleCount == 1 && referenceAlleleIndex == 0 && !canBeReference)) {
            return make_shared<BestAllele<E, A>>(this, sampleIndex, evidenceIndex, -1, -numeric_limits<double>::infinity(), -numeric_limits<double>::infinity());
        }

        auto & sampleValues = valuesBySampleIndex[sampleIndex];
        int bestAlleleIndex = canBeReference || referenceAlleleIndex != 0 ? 0 : 1;

        int secondBestIndex = 0;
        double bestLikelihood = sampleValues[bestAlleleIndex][evidenceIndex];
        double secondBestLikelihood =  -numeric_limits<double>::infinity();

        for(int a = bestAlleleIndex + 1; a < alleleCount; a++){
            if (!canBeReference && referenceAlleleIndex == a) {
                continue;
            }
            double candidateLikelihood = sampleValues[a][evidenceIndex];
            if (candidateLikelihood > bestLikelihood) {
                secondBestIndex = bestAlleleIndex;
                bestAlleleIndex = a;
                secondBestLikelihood = bestLikelihood;
                bestLikelihood = candidateLikelihood;
            } else if (candidateLikelihood > secondBestLikelihood) {
                secondBestIndex = a;
                secondBestLikelihood = candidateLikelihood;
            }
        }

        if (priorities != nullptr && bestLikelihood - secondBestLikelihood < getInformativeThreshold()) {
            double bestPriority = priorities[bestAlleleIndex];
            double secondBestPriority = priorities[secondBestIndex];
            for (int a = 0; a < alleleCount; a++) {
                double candidateLikelihood = sampleValues[a][evidenceIndex];
                if (a == bestAlleleIndex || (!canBeReference && a == referenceAlleleIndex) || bestLikelihood - candidateLikelihood > getInformativeThreshold()) {
                    continue;
                }
                double candidatePriority = priorities[a];

                if (candidatePriority > bestPriority) {
                    secondBestIndex = bestAlleleIndex;
                    bestAlleleIndex = a;
                    secondBestPriority = bestPriority;
                    bestPriority = candidatePriority;
                } else if (candidatePriority > secondBestPriority) {
                    secondBestIndex = a;
                    secondBestPriority = candidatePriority;
                }
            }
        }

        bestLikelihood = sampleValues[bestAlleleIndex][evidenceIndex];
        secondBestLikelihood = secondBestIndex != bestAlleleIndex ? sampleValues[secondBestIndex][evidenceIndex] : -numeric_limits<double>::infinity();
        return make_shared<BestAllele<E, A>>(this, sampleIndex, evidenceIndex, bestAlleleIndex, bestLikelihood, secondBestLikelihood);
    }

    shared_ptr<BestAllele<E, A>> searchBestAllele(int sampleIndex, int evidenceIndex, bool canBeReference){
        return searchBestAllele(sampleIndex, evidenceIndex, canBeReference, nullptr);
    }

    // Does the normalizeLikelihoods job for each piece of evidence.
    void normalizeLikelihoodsPerEvidence(double maximumBestAltLikelihoodDifference, vector<vector<double>>& sampleValues, int sampleIndex, int evidenceIndex){
        //allow the best allele to be the reference because asymmetry leads to strange artifacts like het calls with >90% alt reads
        shared_ptr<BestAllele<E, A>> bestAllele = searchBestAllele(sampleIndex, evidenceIndex, true);

        double worstLikelihoodCap = bestAllele->likelihood + maximumBestAltLikelihoodDifference;

        int alleleCount = alleles->size();

        // Guarantee to be the case by enclosing code.
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][evidenceIndex] < worstLikelihoodCap) {
                sampleValues[a][evidenceIndex] = worstLikelihoodCap;
            }
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

    shared_ptr<vector<shared_ptr<A>>> alleles;


    double maximumLikelihoodOverAllAlleles(int sampleIndex, int evidenceIndex) {
        double result = -numeric_limits<double>::infinity();
        int alleleCount = alleles->size();
        auto & sampleValues = valuesBySampleIndex[sampleIndex];
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][evidenceIndex] > result) {
                result = sampleValues[a][evidenceIndex];
            }
        }
        return result;
    }

public:
    constexpr static double LOG_10_INFORMATIVE_THRESHOLD = 0.2;
    static double NATURAL_LOG_INFORMATIVE_THRESHOLD;

    AlleleLikelihoods(vector<string>& samples, shared_ptr<vector<shared_ptr<A>>>& alleles, map<string, vector<std::shared_ptr<SAMRecord>>>& evidenceBySample) :
         samples(samples), alleles(alleles)
    {
        int sampleCount = samples.size();
        int alleleCount = alleles->size();

        evidenceBySampleIndex.reserve(sampleCount);
        valuesBySampleIndex.reserve(sampleCount);
        referenceAlleleIndex = findReferenceAllele(*alleles);

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

    /**
     * Returns the samples in this evidence-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable. It will not be updated if the collection
     *     allele list changes.
     * </p>
     *
     * @return never {@code null}.
     */
    vector<shared_ptr<A>>& getAlleles(){
        return *alleles;
    }

    /**
     * Adjusts likelihoods so that for each unit of evidence, the best allele likelihood is 0 and caps the minimum likelihood
     * of any allele for each unit of evidence based on the maximum alternative allele likelihood.
     *
     * @param maximumLikelihoodDifferenceCap maximum difference between the best alternative allele likelihood
     *                                           and any other likelihood.
     *
     * @throws IllegalArgumentException if {@code maximumDifferenceWithBestAlternative} is not 0 or less.
     */
    void normalizeLikelihoods(double maximumLikelihoodDifferenceCap){
        assert(!isnan(maximumLikelihoodDifferenceCap) && maximumLikelihoodDifferenceCap < 0.0);

        if(isinf(maximumLikelihoodDifferenceCap))
            return;

        int alleleCount = alleles->size();
        if(alleleCount == 0 || alleleCount == 1)
            return;

        for(int s=0; s<valuesBySampleIndex.size(); s++){
            auto& sampleValues = valuesBySampleIndex[s];
            int evidenceCount = evidenceBySampleIndex[s].size();
            for(int r=0; r<evidenceCount; r++)
            {
                normalizeLikelihoodsPerEvidence(maximumLikelihoodDifferenceCap, sampleValues, s, r);
            }
        }
    }

    /**
   * Removes those read that the best possible likelihood given any allele is just too low.
   *
   * <p>
   *     This is determined by a maximum error per read-base against the best likelihood possible.
   * </p>
   *
   * @param log10MinTrueLikelihood Function that returns the minimum likelihood that the best allele for a unit of evidence must have
   * @throws IllegalStateException is not supported for read-likelihood that do not contain alleles.
   *
   * @throws IllegalArgumentException if {@code maximumErrorPerBase} is negative.
   */   // TODO: validate this method
    void filterPoorlyModeledEvidence(function<double(shared_ptr<SAMRecord>, double)> log10MinTrueLikelihood, double maximumErrorPerBase){
        assert(alleles->size() > 0);
        int numberOfSamples = samples.size();
        for (int s = 0; s < numberOfSamples; s++) {
            auto & sampleEvidence = evidenceBySampleIndex[s];
            vector<int> indexesToRemove;

            int numberOfEvidence = sampleEvidence.size();
            for (int e = 0; e < numberOfEvidence; e++) {
                if (maximumLikelihoodOverAllAlleles(s, e) <  log10MinTrueLikelihood(sampleEvidence[e], maximumErrorPerBase)){
                    indexesToRemove.push_back(e);
                }
            }
            removeSampleEvidence(s, indexesToRemove, alleles->size());
        }
    }

    // ---we use indexes to be removed instead of evidences to be removed
    // Requires that the collection passed iterator can remove elements, and it can be modified.
    void removeSampleEvidence(int sampleIndex, vector<int>& indexesToRemove, int alleleCount) {
        if(indexesToRemove.empty())
            return;

        auto& sampleEvidence = evidenceBySampleIndex[sampleIndex];

        vector<int> IndicesToKeep;
        int numberOfEvidence = sampleEvidence.size();
        int index = 0;
        for(int i=0; i<numberOfEvidence; i++)
        {
            if(index >= indexesToRemove.size())
                break;

            if(indexesToRemove[index] != i)
            {
                IndicesToKeep.push_back(i);
            } else {
                index++;
            }
        }

        // Then we skim out the likelihoods of the removed evidence.
        auto & oldSampleValues = valuesBySampleIndex[sampleIndex];
        vector<vector<double>> newSampleValues(alleleCount, vector<double>(IndicesToKeep.size(), 0.0));

        for (int a = 0; a < alleleCount; a++) {
            for (int r = 0; r < IndicesToKeep.size(); r++) {
                newSampleValues[a][r] = oldSampleValues[a][IndicesToKeep[r]];
            }
        }

        vector<shared_ptr<E>> newSampleEvidence;
        for(int & i : IndicesToKeep)
        {
            newSampleEvidence.template emplace_back(sampleEvidence[i]);
        }

        valuesBySampleIndex[sampleIndex] = newSampleValues;
        evidenceBySampleIndex[sampleIndex] = newSampleEvidence;
    }

    void switchToNaturalLog(){
        assert(!isNaturalLog);
        int sampleCount = samples.size();
        int alleleCount = alleles->size();

        for (int s = 0; s < sampleCount; s++) {
            int evidenceCount = sampleEvidenceCount(s);
            for (int a = 0; a < alleleCount; a++) {
                for (int e = 0; e < evidenceCount; e++) {
                    valuesBySampleIndex[s][a][e] = MathUtils::log10ToLog(valuesBySampleIndex[s][a][e]);
                }
            }
        }
        isNaturalLog = true;
    }

    /**
    * Returns the quantity of evidence that belongs to a sample in the evidence-likelihood collection.
    * @param sampleIndex the query sample index.
    *
    * @return 0 or greater.
    */
    int sampleEvidenceCount(int sampleIndex) {
        assert(sampleIndex >= 0 && sampleIndex < samples.size());
        return evidenceBySampleIndex[sampleIndex].size();
    }

    /**
    * Returns the collection of best allele estimates for the evidence based on the evidence-likelihoods.
    * "Ties" where the ref likelihood is within {@code AlleleLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
    * are broken by the {@code tieBreakingPriority} function.
    *
    * @return never {@code null}, one element per unit fo evidence in the evidence-likelihoods collection.
    */
    shared_ptr<vector<shared_ptr<BestAllele<E, A>>>> bestAllelesBreakingTies(function<double(shared_ptr<A>)> tieBreakingPriority)
    {
        vector<double> priorities;
        if(!alleles->empty())
        {
            for(shared_ptr<A> h : *alleles)
            {
                assert(h != nullptr);
                double tmp = tieBreakingPriority(h);
                priorities.push_back(tmp);

            }

        }

        shared_ptr<vector<shared_ptr<BestAllele<E, A>>>> result = shared_ptr<vector<shared_ptr<BestAllele<E, A>>>>(new vector<shared_ptr<BestAllele<E, A>>>);
        for(int sampleIndex = 0; sampleIndex < samples.size(); sampleIndex++)
        {
            int evidenceCount = evidenceBySampleIndex[sampleIndex].size();

            for(int r=0; r<evidenceCount; r++)
            {
                result->template emplace_back(searchBestAllele(sampleIndex, r, true, priorities.data()));
            }
        }

        return result;
    }

    // TODO: is evidenceIndexBySampleIndex useful ?
    void changeEvidence(shared_ptr<unordered_map<shared_ptr<E>, shared_ptr<E>>> evidenceReplacements)
    {
        int sampleCount = samples.size();
        for(int s = 0; s < sampleCount; s++)
        {
            auto & sampleEvidence = evidenceBySampleIndex[s];
            // Object2IntMap<EVIDENCE> evidenceIndex = evidenceIndexBySampleIndex.get(s);
            int sampleEvidenceCount = sampleEvidence.size();
            for (int r = 0; r < sampleEvidenceCount; r++) {
                shared_ptr<E>& evidence = sampleEvidence[r];
                shared_ptr<E>& replacement = evidenceReplacements->at(evidence);
                if(replacement == nullptr)
                    continue;

                sampleEvidence[r] = replacement;
               /* if (evidenceIndex != null) {
                    evidenceIndex.remove(evidence);
                    evidenceIndex.put(replacement, r);
                }*/
            }
        }
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

    vector<shared_ptr<A>>& alleles(){
        return likelihood->getAlleles();
    }

    AlleleLikelihoods<E, A>& getLikelihoods(){
        assert(likelihood != nullptr);
        return *likelihood;
    }

    int numberOfAlleles(){
        return alleles().size();
    }
};

template <typename E, typename A>
class BestAllele{
public:
    /**
     * Null if there is no possible match (no allele?).
     */
    shared_ptr<A> allele;

    /**
     * The containing sample.
     */
    string sample;

    /**
     * The query evidence.
     */
    shared_ptr<E> evidence;

    /**
     * If allele != null, the indicates the likelihood of the evidence.
     */
    double likelihood;

    AlleleLikelihoods<E, A> * alleleLikelihoods;

    /**
     * Confidence that the evidence actually was generated under that likelihood.
     * This is equal to the difference between this and the second best allele match.
     */
    double confidence;

    BestAllele(AlleleLikelihoods<E, A> * likelihoods, int sampleIndex, int evidenceIndex, int bestAlleleIndex,
               double likelihood, double secondBestLikelihood){
        allele = bestAlleleIndex == -1 ? nullptr : likelihoods->alleles->operator[](bestAlleleIndex);
        this->likelihood = likelihood;
        sample = likelihoods->samples[sampleIndex];
        evidence = likelihoods->evidenceBySampleIndex[sampleIndex][evidenceIndex];
        confidence = likelihood == secondBestLikelihood ? 0 : likelihood - secondBestLikelihood;
        alleleLikelihoods = likelihoods;
    }

    bool isInformative() {
        return confidence > AlleleLikelihoods<E, A>::LOG_10_INFORMATIVE_THRESHOLD;
    }
};

#endif //MUTECT2CPP_MASTER_ALLELELIKELIHOODS_H
