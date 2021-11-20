//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
#define MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H

#include "AssemblyResult.h"
#include "haplotype/Haplotype.h"

class ReadThreadingAssembler {
private:
    bool recoverDanglingBranches = true;
    int pruneFactor;
    int minDanglingBranchLength = 0;
    bool recoverAllDanglingBranches = false;
    AssemblyResult* getAssemblyResult(Haplotype* refHaplotype, int kmerSize, ReadThreadingGraph* rtgraph);
    AssemblyResult* cleanupSeqGraph(SeqGraph* seqGraph);
};


#endif //MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
