//
// Created by 梦想家xixi on 2021/10/18.
//

#ifndef MUTECT2CPP_MASTER_READTHREADINGGRAPH_H
#define MUTECT2CPP_MASTER_READTHREADINGGRAPH_H

#include "SAMRecord.h"
#include "MultiDeBruijnVertex.h"
#include "MultiSampleEdge.h"
#include <map>
#include <string>
#include <vector>
#include "Kmer.h"
#include "BaseGraph/DirectedSpecifics.h"

typedef struct{
    std::string name;
    uint8_t* sequence;
    int start;
    int stop;
    int count;
    bool isRef;
}SequenceForKmers;

class ReadThreadingGraph : public DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>{
protected:
    int kmerSize;

private:
    int numPruningSamples;

    Kmer refSource;

    ArraySet<Kmer> nonUniqueKmers;

    std::map<Kmer, MultiDeBruijnVertex*> uniqueKmers;

    const uint8_t minBaseQualityToUseInAssembly;

    bool startThreadingOnlyAtExistingVertex;

    bool increaseCountsThroughBranches = false;

    static const bool INCREASE_COUNTS_BACKWARDS = true;
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
     ArraySet<Kmer> determineKmerSizeAndNonUniques(int minKmerSize, int maxKmerSize);

    /**
    * Create a new vertex for kmer.  Add it to the uniqueKmers map if appropriate.
    *
    * kmer must not have a entry in unique kmers, or an error will be thrown
    *
    * @param kmer the kmer we want to create a vertex for
    * @return the non-null created vertex
    */
     MultiDeBruijnVertex* createVertex(Kmer & kmer);

    /**
   * Workhorse routine of the assembler.  Given a sequence whose last vertex is anchored in the graph, extend
   * the graph one bp according to the bases in sequence.
   *
   * @param prevVertex a non-null vertex where sequence was last anchored in the graph
   * @param sequence the sequence we're threading through the graph
   * @param kmerStart the start of the current kmer in graph we'd like to add
   * @param count the number of observations of this kmer in graph (can be > 1 for GGA)
   * @param isRef is this the reference sequence?
   * @return a non-null vertex connecting prevVertex to in the graph based on sequence
   */
     MultiDeBruijnVertex* extendChainByOne(MultiDeBruijnVertex* prevVertex, uint8_t * sequence, int kmerStart, int count, bool isRef);

     void threadSequence(SequenceForKmers & sequenceForKmers);

    /**
    * Find vertex and its position in seqForKmers where we should start assembling seqForKmers
    *
    * @param seqForKmers the sequence we want to thread into the graph
    * @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
    */
     int findStart(SequenceForKmers seqForKmers);

     bool getUniqueKmerVertex(Kmer & kmer, bool allowRefSource);

     MultiDeBruijnVertex* getOrCreateKmerVertex(uint8_t * sequence, int start);

     void increaseCountsInMatchedKmers(SequenceForKmers & seqForKmers, MultiDeBruijnVertex* vertex, uint8_t* originalKmer, int offset);

public:
    ReadThreadingGraph(uint8_t minBaseQualityToUseInAssembly, int kmerSize, bool alreadyBuilt, Kmer ref, int numPruningSamples) : minBaseQualityToUseInAssembly(minBaseQualityToUseInAssembly), kmerSize(kmerSize), alreadyBuilt(
            false), refSource(ref), numPruningSamples(numPruningSamples){}

    static std::vector<Kmer> determineNonUniqueKmers(SequenceForKmers &sequenceForKmers, int kmerSize);

    /**
     * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
     * been added to the graph.
     */
    void buildGraphIfNecessary();

    /**
     * Changes the threading start location policy.
     *
     * @param value  {@code true} if threading will start only at existing vertices in the graph, {@code false} if
     *  it can start at any unique kmer.
     */
     void setThreadingStartOnlyAtExistingVertex(bool value);

    /**
     * Add a read to the sequence graph.  Finds maximal consecutive runs of bases with sufficient quality
     * and applies {@see addSequence} to these subreads if they are longer than the kmer size.
     *
     * @param read a non-null read
     */
     void addRead(SAMRecord read);

     bool removeVertex(MultiDeBruijnVertex* V);

     void setPending();


    void removeSingletonOrphanVertices();
};


#endif //MUTECT2CPP_MASTER_READTHREADINGGRAPH_H
