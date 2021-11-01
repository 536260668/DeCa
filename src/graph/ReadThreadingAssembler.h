//
// Created by 梦想家xixi on 2021/10/18.
//

#ifndef MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
#define MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H

#include "SAMRecord.h"
#include "MultiDeBruijnVertex.h"
#include <map>
#include <set>
#include <string>
#include <vector>
#include "Kmer.h"

typedef struct{
    std::string name;
    uint8_t* sequence;
    int start;
    int stop;
    int count;
    bool isRef;
}SequenceForKmers;

class ReadThreadingAssembler {
protected:
    /**
     * Add a read to the sequence graph.  Finds maximal consecutive runs of bases with sufficient quality
     * and applies {@see addSequence} to these subreads if they are longer than the kmer size.
     *
     * @param read a non-null read
     */
    void addRead(SAMRecord read);
    int kmerSize;

private:
    std::set<MultiDeBruijnVertex> unmodifiableVertexSet;

    std::set<Kmer> nonUniqueKmers;

    std::map<Kmer, MultiDeBruijnVertex> uniqueKmers;

    const uint8_t minBaseQualityToUseInAssembly;
    /**
     * Determines whether a base can safely be used for assembly.
     * Currently disallows Ns and/or those with low quality
     *
     * @param base  the base under consideration
     * @param qual  the quality of that base
     * @return true if the base can be used for assembly, false otherwise
     */
    bool baseIsUsableForAssembly(uint8_t base, uint8_t qual) const;

    bool alreadyBuilt;

    std::map<std::string, std::vector<SequenceForKmers>> pending;

    void addSequence(std::string seqName, std::string& sampleName, const uint8_t* sequence, int start, int stop, int count, bool isRef);

    /**
     * Get the collection of all sequences for kmers across all samples in no particular order
     * @return non-null Collection
     */
    std::vector<SequenceForKmers> getAllPendingSequences();

    /**
     * Compute the smallest kmer size >= minKmerSize and <= maxKmerSize that has no non-unique kmers
     * among all sequences added to the current graph.  Will always return a result for maxKmerSize if
     * all smaller kmers had non-unique kmers.
     *
     * @param minKmerSize the minimum kmer size to consider when constructing the graph
     * @param maxKmerSize the maximum kmer size to consider
     * @return a non-null NonUniqueResult
     */
     std::set<Kmer> determineKmerSizeAndNonUniques(int minKmerSize, int maxKmerSize);

    /**
    * Create a new vertex for kmer.  Add it to the uniqueKmers map if appropriate.
    *
    * kmer must not have a entry in unique kmers, or an error will be thrown
    *
    * @param kmer the kmer we want to create a vertex for
    * @return the non-null created vertex
    */
     MultiDeBruijnVertex createVertex(Kmer kmer);

public:
    ReadThreadingAssembler(uint8_t minBaseQualityToUseInAssembly, int kmerSize, bool alreadyBuilt) : minBaseQualityToUseInAssembly(minBaseQualityToUseInAssembly), kmerSize(kmerSize), alreadyBuilt(
            false) {}

    static std::vector<Kmer> determineNonUniqueKmers(SequenceForKmers &sequenceForKmers, int kmerSize);

    void test(Kmer kmer);

};


#endif //MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
