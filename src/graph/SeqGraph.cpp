//
// Created by 梦想家xixi on 2021/11/15.
//

#include "SeqGraph.h"
#include <list>
#include <climits>
#include <cstring>
#include "utils/MergeCommonSuffices.h"
#include "utils/MergeDiamonds.h"
#include "utils/MergeTails.h"
#include "utils/SharedSequenceMerger.h"
#include "utils/SplitCommonSuffices.h"
#include "utils/GraphUtils.h"


std::shared_ptr<BaseEdge> SeqGraph::createEdge(std::shared_ptr<SeqVertex> sourceVertex, std::shared_ptr<SeqVertex> targetVertrx) {
    return std::shared_ptr<BaseEdge>(new BaseEdge(false, 1));
}

bool SeqGraph::zipLinearChains() {
    std::list<std::shared_ptr<SeqVertex>> zipStarts;
    for(const std::shared_ptr<SeqVertex>& source : DirectedSpecifics<SeqVertex, BaseEdge>::getVertexSet()) {
        if(isLinearChainStart(source)) {
            zipStarts.emplace_back(source);
        }
    }

    if(zipStarts.empty())
        return false;

    bool mergedOne = false;
    for(const std::shared_ptr<SeqVertex>& zipStart : zipStarts) {
        std::list<std::shared_ptr<SeqVertex>> linearChain = traceLinearChain(zipStart);

        mergedOne |= mergeLinearChain(linearChain);
    }
    return mergedOne;
}

bool SeqGraph::isLinearChainStart(std::shared_ptr<SeqVertex>source) {
    return DirectedSpecifics<SeqVertex, BaseEdge>::outDegreeOf(source) == 1
    && (DirectedSpecifics<SeqVertex, BaseEdge>::inDegreeOf(source) != 1 ||
            DirectedSpecifics<SeqVertex, BaseEdge>::outDegreeOf(*(incomingVerticesOf(source).begin())) > 1);
}

std::list<std::shared_ptr<SeqVertex>> SeqGraph::traceLinearChain(std::shared_ptr<SeqVertex>zipStart) {
    std::list<std::shared_ptr<SeqVertex>> linearChain;
    linearChain.emplace_back(zipStart);

    bool lastIsRef = isReferenceNode(zipStart);
    std::shared_ptr<SeqVertex> last = zipStart;
    while(true) {
        if (DirectedSpecifics<SeqVertex, BaseEdge>::outDegreeOf(last) != 1) {
            break;
        }

        std::shared_ptr<SeqVertex> target = getEdgeTarget(outgoingEdgeOf(last));

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

bool SeqGraph::mergeLinearChain(std::list<std::shared_ptr<SeqVertex>> &linearChain) {
    Mutect2Utils::validateArg(!linearChain.empty(), "BUG: cannot have linear chain with 0 elements");

    std::shared_ptr<SeqVertex> first = linearChain.front();
    std::shared_ptr<SeqVertex> last = linearChain.back();

    if(first == last)
        return false;

    std::shared_ptr<SeqVertex> addedVertex = mergeLinearChainVertices(linearChain);
    DirectedSpecifics<SeqVertex, BaseEdge>::addVertex(addedVertex);

    for(const std::shared_ptr<BaseEdge>& edge : DirectedSpecifics<SeqVertex, BaseEdge>::outgoingEdgesOf(last)) {
        addEdge(addedVertex, getEdgeTarget(edge), std::shared_ptr<BaseEdge>(new BaseEdge(edge->getIsRef(), edge->getMultiplicity())));
    }

    for(const std::shared_ptr<BaseEdge>& edge : DirectedSpecifics<SeqVertex, BaseEdge>::incomingEdgesOf(first)) {
        addEdge(getEdgeSource(edge), addedVertex, std::shared_ptr<BaseEdge>(new BaseEdge(edge->getIsRef(), edge->getMultiplicity())));
    }
    DirectedSpecifics<SeqVertex, BaseEdge>::removeAllVertices(linearChain);
    return true;
}

std::shared_ptr<SeqVertex> SeqGraph::mergeLinearChainVertices(std::list<std::shared_ptr<SeqVertex>> &vertices) {
    int length = 500;
    int start = 0;
    std::shared_ptr<uint8_t[]> tmp(new uint8_t[length]);
    for(std::shared_ptr<SeqVertex> v : vertices) {
        int seqLength = v->getLength();
        std::shared_ptr<uint8_t[]> seq = v->getSequence();
        while(start + seqLength >= length) {
            length *= 2;
            std::shared_ptr<uint8_t[]> newtmp(new uint8_t[length]);
            memcpy(newtmp.get(), tmp.get(), start);
            tmp = newtmp;
        }
        memcpy(tmp.get()+start, seq.get(), seqLength);
        start += seqLength;
    }
    return std::shared_ptr<SeqVertex>(new SeqVertex(tmp, start));
}

void SeqGraph::simplifyGraph() {
    simplifyGraph(INT_MAX);
}

void SeqGraph::simplifyGraph(int maxCycles) {
    zipLinearChains();
    std::shared_ptr<SeqGraph> prevGraph = nullptr;
    for(int i = 0; i < maxCycles; i++) {
        if (i > MAX_REASONABLE_SIMPLIFICATION_CYCLES) {
            throw std::invalid_argument("Infinite loop detected in simplification routines for kmer graph");
        }
        bool didSomeWork = simplifyGraphOnce(i);
        if(! didSomeWork) {
            break;
        }
        if(i > 5) {
            if(prevGraph != nullptr && GraphUtils::graphEquals(prevGraph.get(), this))
                break;
        }
        prevGraph = std::shared_ptr<SeqGraph>(new SeqGraph(*this));
    }
}

bool SeqGraph::simplifyGraphOnce(int iteration) {
    bool didSomeWork = false;
    std::shared_ptr<SeqGraph> graph(new SeqGraph(*this));
    didSomeWork |= MergeDiamonds(graph).transformUntilComplete();
    didSomeWork |= MergeTails(graph).transformUntilComplete();
    didSomeWork |= SplitCommonSuffices(graph).transformUntilComplete();
    didSomeWork |= MergeCommonSuffices(graph).transformUntilComplete();
    didSomeWork |= zipLinearChains();
    return didSomeWork;
}

SeqGraph::SeqGraph(SeqGraph &seqGraph) : kmerSize(seqGraph.kmerSize), DirectedSpecifics<SeqVertex, BaseEdge>(){
    vertexMapDirected = seqGraph.vertexMapDirected;
    edgeMap = seqGraph.edgeMap;
}

std::shared_ptr<SeqGraph> SeqGraph::clone() {
    std::shared_ptr<SeqGraph> ret(new SeqGraph(*this));
    return ret;
}

