//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
#define MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H

#include "AssemblyResult.h"
#include "haplotype/Haplotype.h"
#include "SimpleInterval.h"

class ReadThreadingAssembler {
private:
    bool recoverDanglingBranches = true;
    int pruneFactor;
    int minDanglingBranchLength = 0;
    bool recoverAllDanglingBranches = false;
    bool removePathsNotConnectedToRef = true;
    bool justReturnRawGraph = false;
    AssemblyResult* getAssemblyResult(Haplotype* refHaplotype, int kmerSize, ReadThreadingGraph* rtgraph);
    AssemblyResult* cleanupSeqGraph(SeqGraph* seqGraph);
    std::vector<Haplotype*> findBestPaths(std::list<SeqGraph*> graph, Haplotype* refHaplotype, SimpleInterval* refLoc, SimpleInterval* activeRegionWindow,
                                          std::map<SeqGraph*, AssemblyResult*> assemblyResultByGraph);
};


#endif //MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
