//
// Created by 梦想家xixi on 2021/11/20.
//

#include "SharedSequenceMerger.h"

bool SharedSequenceMerger::canMerge(SeqGraph *graph, SeqVertex *v, ArraySet<SeqVertex *> incomingVertices) {
    if(incomingVertices.empty()) {
        return false;
    }

    SeqVertex* first = *incomingVertices.begin();
    for(SeqVertex* prev : incomingVertices) {
        if(! prev->seqEquals(first)) {
            return false;
        }
        ArraySet<SeqVertex*> prevOuts = graph->outgoingVerticesOf(prev);
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

bool SharedSequenceMerger::merge(SeqGraph *graph, SeqVertex *v) {
    Mutect2Utils::validateArg(graph, "graph cannot be null");
    ArraySet<SeqVertex*> allVertex = graph->getVertexSet();
    Mutect2Utils::validateArg(allVertex.find(v) != allVertex.end(), "graph doesn't contain vertex");
    ArraySet<SeqVertex*> prevs = graph->incomingVerticesOf(v);
    if(!canMerge(graph, v, prevs)) {
        return false;
    } else {
        std::list<BaseEdge*> edgesToRemove;
        uint8_t * prevSeq = (*prevs.begin())->getSequence();
        int prevSeqLength = (*prevs.begin())->getLength();
        uint8_t * vSeq = v->getSequence();
        int vSeqLength = v->getLength();
        int tmpLength = prevSeqLength + vSeqLength;
        uint8_t * tmp = new uint8_t[tmpLength];
        memcpy(tmp, prevSeq, prevSeqLength);
        memcpy(tmp+prevSeqLength, vSeq, vSeqLength);
        SeqVertex* newV = new SeqVertex(tmp, tmpLength);
        graph->addVertex(newV);
        for(SeqVertex* prev : prevs) {
            for(BaseEdge* prevIn : graph->incomingEdgesOf(prev)) {
                graph->addEdge(graph->getEdgeSource(prevIn), newV, new BaseEdge(prevIn->getIsRef(), prevIn->getMultiplicity()));
                edgesToRemove.emplace_back(prevIn);
            }
        }
        for(BaseEdge* e : graph->outgoingEdgesOf(v)) {
            graph->addEdge(newV, graph->getEdgeTarget(e), new BaseEdge(e->getIsRef(), e->getMultiplicity()));
        }
        graph->removeAllVertices(prevs.getArraySet());
        graph->removeVertex(v);
        graph->removeAllEdges(edgesToRemove);
        return true;
    }
}
