//
// Created by 梦想家xixi on 2021/11/15.
//

#include "SeqGraph.h"
#include <list>
#include <climits>
#include "utils/MergeCommonSuffices.h"
#include "utils/MergeDiamonds.h"
#include "utils/MergeTails.h"
#include "utils/SharedSequenceMerger.h"
#include "utils/SplitCommonSuffices.h"
#include "utils/GraphUtils.h"


BaseEdge *SeqGraph::createEdge(SeqVertex *sourceVertex, SeqVertex *targetVertrx) {
    return new BaseEdge(false, 1);
}

bool SeqGraph::zipLinearChains() {
    std::list<SeqVertex*> zipStarts;
    for(SeqVertex* source : DirectedSpecifics<SeqVertex, BaseEdge>::getVertexSet()) {
        if(isLinearChainStart(source)) {
            zipStarts.emplace_back(source);
        }
    }

    if(zipStarts.empty())
        return false;

    bool mergedOne = false;
    for(SeqVertex* zipStart : zipStarts) {
        std::list<SeqVertex*> linearChain = traceLinearChain(zipStart);

        mergedOne |= mergeLinearChain(linearChain);
    }
    return mergedOne;
}

bool SeqGraph::isLinearChainStart(SeqVertex *source) {
    return DirectedSpecifics<SeqVertex, BaseEdge>::outDegreeOf(source) == 1
    && (DirectedSpecifics<SeqVertex, BaseEdge>::inDegreeOf(source) != 1 ||
            DirectedSpecifics<SeqVertex, BaseEdge>::outDegreeOf(*(incomingVerticesOf(source).begin())) > 1);
}

std::list<SeqVertex *> SeqGraph::traceLinearChain(SeqVertex *zipStart) {
    std::list<SeqVertex *> linearChain;
    linearChain.emplace_back(zipStart);

    bool lastIsRef = isReferenceNode(zipStart);
    SeqVertex* last = zipStart;
    while(true) {
        if (DirectedSpecifics<SeqVertex, BaseEdge>::outDegreeOf(last) != 1) {
            break;
        }

        SeqVertex* target = getEdgeTarget(outgoingEdgeOf(last));

        if(DirectedSpecifics<SeqVertex, BaseEdge>::inDegreeOf(target) != 1 || last == target) {
            break;
        }

        bool targetIsRef = isReferenceNode(target);
        if (lastIsRef != targetIsRef) {
            break;
        }
        linearChain.emplace_back(target);
        last = target;
        lastIsRef = targetIsRef;
    }
    return linearChain;
}

bool SeqGraph::mergeLinearChain(std::list<SeqVertex *> &linearChain) {
    Mutect2Utils::validateArg(!linearChain.empty(), "BUG: cannot have linear chain with 0 elements");

    SeqVertex* first = linearChain.front();
    SeqVertex* last = linearChain.back();

    if(first == last)
        return false;

    SeqVertex* addedVertex = mergeLinearChainVertices(linearChain);
    DirectedSpecifics<SeqVertex, BaseEdge>::addVertex(addedVertex);

    for(BaseEdge* edge : DirectedSpecifics<SeqVertex, BaseEdge>::outgoingEdgesOf(last)) {
        addEdge(addedVertex, getEdgeTarget(edge), new BaseEdge(edge->getIsRef(), edge->getMultiplicity()));
    }

    for(BaseEdge* edge : DirectedSpecifics<SeqVertex, BaseEdge>::incomingEdgesOf(first)) {
        addEdge(getEdgeSource(edge), addedVertex, new BaseEdge(edge->getIsRef(), edge->getMultiplicity()));
    }
    DirectedSpecifics<SeqVertex, BaseEdge>::removeAllVertices(linearChain);
    return true;
}

SeqVertex *SeqGraph::mergeLinearChainVertices(std::list<SeqVertex *> &vertices) {
    int length = 500;
    int start = 0;
    uint8_t * tmp = new uint8_t[length];
    for(SeqVertex* v : vertices) {
        int seqLength = v->getLength();
        uint8_t * seq = v->getSequence();
        while(start + seqLength >= length) {
            length *= 2;
            uint8_t * newtmp = new uint8_t[length];
            memcpy(newtmp, tmp, start);
            delete[] tmp;
            tmp = newtmp;
        }
        memcpy(tmp+start, seq, seqLength);
        start += seqLength;
    }
    return new SeqVertex(tmp, start);
}

void SeqGraph::simplifyGraph() {
    simplifyGraph(INT_MAX);
}

void SeqGraph::simplifyGraph(int maxCycles) {
    zipLinearChains();
    SeqGraph* prevGraph = nullptr;
    for(int i = 0; i < maxCycles; i++) {
        if (i > MAX_REASONABLE_SIMPLIFICATION_CYCLES) {
            throw std::invalid_argument("Infinite loop detected in simplification routines for kmer graph");
        }
        bool didSomeWork = simplifyGraphOnce(i);
        if(! didSomeWork) {
            break;
        }
        if(i > 5) {
            if(prevGraph != nullptr && GraphUtils::graphEquals(prevGraph, this))
                break;
        }
        prevGraph = new SeqGraph(*this);
    }
}

bool SeqGraph::simplifyGraphOnce(int iteration) {
    bool didSomeWork = false;
    didSomeWork |= MergeDiamonds(this).transformUntilComplete();
    didSomeWork |= MergeTails(this).transformUntilComplete();
    didSomeWork |= SplitCommonSuffices(this).transformUntilComplete();
    didSomeWork |= MergeCommonSuffices(this).transformUntilComplete();
    didSomeWork |= zipLinearChains();
    return didSomeWork;
}

SeqGraph::SeqGraph(SeqGraph &seqGraph) : kmerSize(seqGraph.kmerSize), DirectedSpecifics<SeqVertex, BaseEdge>(){
    vertexMapDirected = seqGraph.vertexMapDirected;
    edgeMap = seqGraph.edgeMap;
}

SeqGraph *SeqGraph::clone() {
    SeqGraph* ret = new SeqGraph(*this);
    return ret;
}

