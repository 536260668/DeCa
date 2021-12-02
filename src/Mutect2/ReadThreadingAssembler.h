//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
#define MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H

#include "AssemblyResult.h"
#include "haplotype/Haplotype.h"
#include "SimpleInterval.h"
#include "AssemblyResultSet.h"
#include "tools/haplotypecaller/ReadErrorCorrector.h"
#include "path/ChainPruner.h"

class ReadThreadingAssembler {
private:
    bool recoverDanglingBranches = true;
    int pruneFactor;
    int numPruningSamples;
    int minDanglingBranchLength = 0;
    int numBestHaplotypesPerGraph;
    bool recoverAllDanglingBranches = false;
    bool removePathsNotConnectedToRef = true;
    bool justReturnRawGraph = false;
    bool debugGraphTransformations = false;
    ChainPruner<MultiDeBruijnVertex, MultiSampleEdge>* chainPruner;
    static const uint8_t DEFAULT_MIN_BASE_QUALITY_TO_USE = 10;
    AssemblyResult* getAssemblyResult(Haplotype* refHaplotype, int kmerSize, ReadThreadingGraph* rtgraph);
    AssemblyResult* cleanupSeqGraph(SeqGraph* seqGraph);
    std::vector<Haplotype*> findBestPaths(const std::list<SeqGraph*>& graph, Haplotype* refHaplotype, SimpleInterval* refLoc, SimpleInterval* activeRegionWindow,
                                          const std::map<SeqGraph*, AssemblyResult*>& assemblyResultByGraph, AssemblyResultSet* assemblyResultSet) const;
    AssemblyResult* createGraph(std::vector<SAMRecord*> reads, Haplotype* refHaplotype, int kmerSize, bool allowLowComplexityGraphs, bool allowNonUniqueKmersInRef);

public:
    AssemblyResultSet* runLocalAssembly(AssemblyRegion * assemblyRegion, Haplotype* refHaplotype, uint8_t* fullReferenceWithPadding, int refLength, SimpleInterval* refLoc, ReadErrorCorrector* readErrorCorrector);
    std::vector<AssemblyResult*> assemble(std::vector<SAMRecord> & reads, Haplotype* refHaplotype);

protected:
    uint8_t minBaseQualityToUseInAssembly = DEFAULT_MIN_BASE_QUALITY_TO_USE;
};


#endif //MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
