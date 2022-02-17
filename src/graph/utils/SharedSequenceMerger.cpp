//
// Created by 梦想家xixi on 2021/11/20.
//

#include "SharedSequenceMerger.h"

bool SharedSequenceMerger::canMerge(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> v, ArraySet<std::shared_ptr<SeqVertex>> incomingVertices) {
    if(incomingVertices.empty()) {
        return false;
    }

    std::shared_ptr<SeqVertex> first = *incomingVertices.begin();
    for(std::shared_ptr<SeqVertex> prev : incomingVertices) {
        if(! prev->seqEquals(first)) {
            return false;
        }
        ArraySet<std::shared_ptr<SeqVertex>> prevOuts = graph->outgoingVerticesOf(prev);
        if(prevOuts.size() != 1){
            return false;
        }
        if(*prevOuts.begin() != v) {
            return false;
        }
        if(graph->inDegreeOf(prev) == 0) {
            return false;
        }
    }
    return true;
}

bool SharedSequenceMerger::merge(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> v) {
    Mutect2Utils::validateArg(graph.get(), "graph cannot be null");
    ArraySet<std::shared_ptr<SeqVertex>> allVertex = graph->getVertexSet();
    Mutect2Utils::validateArg(allVertex.find(v) != allVertex.end(), "graph doesn't contain vertex");
    ArraySet<std::shared_ptr<SeqVertex>> prevs = graph->incomingVerticesOf(v);
    if(!canMerge(graph, v, prevs)) {
        return false;
    } else {
        std::list<std::shared_ptr<BaseEdge>> edgesToRemove;
        std::shared_ptr<uint8_t[]> prevSeq = (*prevs.begin())->getSequence();
        int prevSeqLength = (*prevs.begin())->getLength();
        std::shared_ptr<uint8_t[]> vSeq = v->getSequence();
        int vSeqLength = v->getLength();
        int tmpLength = prevSeqLength + vSeqLength;
        std::shared_ptr<uint8_t[]> tmp(new uint8_t[tmpLength]);
        memcpy(tmp.get(), prevSeq.get(), prevSeqLength);
        memcpy(tmp.get()+prevSeqLength, vSeq.get(), vSeqLength);
        std::shared_ptr<SeqVertex> newV(new SeqVertex(tmp, tmpLength));
        graph->addVertex(newV);
        for(std::shared_ptr<SeqVertex> prev : prevs) {
            for(std::shared_ptr<BaseEdge> prevIn : graph->incomingEdgesOf(prev)) {
                graph->addEdge(graph->getEdgeSource(prevIn), newV, std::shared_ptr<BaseEdge>(new BaseEdge(prevIn->getIsRef(), prevIn->getMultiplicity())));
                edgesToRemove.emplace_back(prevIn);
            }
        }
        for(std::shared_ptr<BaseEdge> e : graph->outgoingEdgesOf(v)) {
            graph->addEdge(newV, graph->getEdgeTarget(e), std::shared_ptr<BaseEdge>(new BaseEdge(e->getIsRef(), e->getMultiplicity())));
        }
        graph->removeAllVertices(prevs.getArraySet());
        graph->removeVertex(v);
        graph->removeAllEdges(edgesToRemove);
        return true;
    }
}
