//
// Created by 梦想家xixi on 2021/11/15.
//

#include "ReadThreadingAssembler.h"

AssemblyResult *
ReadThreadingAssembler::getAssemblyResult(Haplotype *refHaplotype, int kmerSize, ReadThreadingGraph *rtgraph) {
    if(recoverDanglingBranches){
        rtgraph->recoverDanglingTails(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
        rtgraph->recoverDanglingHeads(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
    }
}

AssemblyResult *ReadThreadingAssembler::cleanupSeqGraph(SeqGraph *seqGraph) {
    seqGraph->zipLinearChains();

}
