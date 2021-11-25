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

    if(removePathsNotConnectedToRef){
        rtgraph->removePathsNotConnectedToRef();
    }

    SeqGraph* initialSeqGraph = rtgraph->toSequenceGraph();
    if(justReturnRawGraph) {
        return new AssemblyResult(ASSEMBLED_SOME_VARIATION, initialSeqGraph, nullptr);
    }
    initialSeqGraph->cleanNonRefPaths();

    AssemblyResult* cleaned = cleanupSeqGraph(initialSeqGraph);
    Status status = cleaned->getStatus();
    AssemblyResult* ret = new AssemblyResult(status, cleaned->getGraph(), rtgraph);
    delete cleaned;
    return ret;
}

AssemblyResult *ReadThreadingAssembler::cleanupSeqGraph(SeqGraph *seqGraph) {
    seqGraph->zipLinearChains();
    seqGraph->removeSingletonOrphanVertices();
    seqGraph->removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
    seqGraph->simplifyGraph();
    if(seqGraph->getReferenceSinkVertex() == nullptr || seqGraph->getReferenceSourceVertex() == nullptr) {
        return new AssemblyResult(JUST_ASSEMBLED_REFERENCE, seqGraph, nullptr);
    }
    seqGraph->removePathsNotConnectedToRef();
    seqGraph->simplifyGraph();
    if(seqGraph->getVertexSet().size() == 1) {
        SeqVertex* complete = *(seqGraph->getVertexSet().begin());
        SeqVertex* dummy = new SeqVertex(nullptr, 0);
        seqGraph->addVertex(dummy);
        seqGraph->addEdge(complete, dummy, new BaseEdge(true, 0));
    }
    return new AssemblyResult(ASSEMBLED_SOME_VARIATION, seqGraph, nullptr);
}
