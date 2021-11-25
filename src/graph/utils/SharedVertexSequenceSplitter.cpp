//
// Created by 梦想家xixi on 2021/11/19.
//

#include "SharedVertexSequenceSplitter.h"
#include "GraphUtils.h"

std::pair<SeqVertex *, SeqVertex *>*
SharedVertexSequenceSplitter::commonPrefixAndSuffixOfVertices(ArraySet<SeqVertex *> middleVertices){
    std::list<std::pair<uint8_t *, int>> kmers;
    int min = INT32_MAX;
    for(SeqVertex* v : middleVertices) {
        std::pair<uint8_t*, int> tmp;
        tmp.first = v->getSequence();
        tmp.second = v->getLength();
        kmers.emplace_back(tmp);
        min = std::min(min, tmp.second);
    }

    int prefixLen = GraphUtils::commonMaximumPrefixLength(kmers);
    int suffixLen = GraphUtils::commonMaximumSuffixLength(kmers, min - prefixLen);

    uint8_t * kmer = kmers.begin()->first;
    int length = kmers.begin()->second;
    int prefixLength;
    uint8_t * prefix = Mutect2Utils::copyOfRange(kmer, length, 0, prefixLen, prefixLength);
    int suffixLength;
    uint8_t * suffix = Mutect2Utils::copyOfRange(kmer, length, length - suffixLen, length, suffixLength);
    return  new std::pair<SeqVertex *, SeqVertex *>(new SeqVertex(prefix, prefixLength), new SeqVertex(suffix, suffixLength));
}

SharedVertexSequenceSplitter::SharedVertexSequenceSplitter(SeqGraph *graph, ArraySet<SeqVertex *> toSplitsArg) : outer(graph), toSplits(toSplitsArg){
    Mutect2Utils::validateArg(graph, "graph cannot be null");
    Mutect2Utils::validateArg(toSplitsArg.size() > 1, "Can only split at least 2 vertices");
    for(SeqVertex* v : toSplitsArg) {
        ArraySet<SeqVertex*> allVertex = graph->getVertexSet();
        if(allVertex.find(v) == allVertex.end())
            throw std::invalid_argument("graph doesn't contain all of the vertices to split");
    }
    std::pair<SeqVertex*, SeqVertex*>* prefixAndSuffix = commonPrefixAndSuffixOfVertices(toSplits);
    prefixV = prefixAndSuffix->first;
    suffixV = prefixAndSuffix->second;
    delete prefixAndSuffix;
}

bool SharedVertexSequenceSplitter::meetsMinMergableSequenceForEitherPrefixOrSuffix(int minCommonSequence) {
    return meetsMinMergableSequenceForPrefix(minCommonSequence) || meetsMinMergableSequenceForSuffix(minCommonSequence);
}

bool SharedVertexSequenceSplitter::meetsMinMergableSequenceForPrefix(const int minCommonSequence) {
    return getPrefixV()->getLength() >= minCommonSequence;
}

bool SharedVertexSequenceSplitter::meetsMinMergableSequenceForSuffix(int minCommonSequence) {
    return getSuffixV()->getLength() >= minCommonSequence;
}

bool SharedVertexSequenceSplitter::splitAndUpdate(SeqVertex *top, SeqVertex *bottom) {
    split();
    updateGraph(top, bottom);
    return true;
}

void SharedVertexSequenceSplitter::split() {
    splitGraph = new SeqGraph(outer->getKmerSize());
    splitGraph->addVertex(getPrefixV());
    splitGraph->addVertex(getSuffixV());
    for(SeqVertex* mid : toSplits) {
        BaseEdge* toMid = processEdgeToRemove(mid, outer->incomingEdgeOf(mid));
        BaseEdge* fromMid = processEdgeToRemove(mid, outer->outgoingEdgeOf(mid));

        SeqVertex* remaining = mid->withoutPrefixAndSuffix(getPrefixV()->getSequence(), getPrefixV()->getLength(), getSuffixV()->getSequence(), getSuffixV()->getLength());
        if(remaining != nullptr) {
            splitGraph->addVertex(remaining);
            getNewMiddles().emplace_back(remaining);
            splitGraph->addEdge(getPrefixV(), remaining, toMid);
            splitGraph->addEdge(remaining, getSuffixV(), fromMid);
        } else {
            BaseEdge* tmp = new BaseEdge(toMid->getIsRef(), toMid->getMultiplicity());
            tmp->add(*fromMid);
            splitGraph->addOrUpdateEdge(getPrefixV(), getSuffixV(), tmp);
        }
    }
}

BaseEdge *SharedVertexSequenceSplitter::processEdgeToRemove(SeqVertex *v, BaseEdge *e) {
    if (e == nullptr) {
        return new BaseEdge(outer->isReferenceNode(v), 0);
    } else {
        edgesToRemove.emplace_back(e);
        return new BaseEdge(e->getIsRef(), e->getMultiplicity());
    }
}

void SharedVertexSequenceSplitter::updateGraph(SeqVertex *top, SeqVertex *bot) {
    for(SeqVertex* v : toSplits) {
        ArraySet<SeqVertex*> allVertex = outer->getVertexSet();
        if(allVertex.find(v) == allVertex.end())
            throw std::invalid_argument("graph doesn't contain all of the vertices to split");
    }
    Mutect2Utils::validateArg(top != nullptr || bot != nullptr, "Cannot update graph without at least one top or bot vertex, but both were null");
    Mutect2Utils::validateArg(top == nullptr || outer->containsVertex(top), "top not found in graph");
    Mutect2Utils::validateArg(bot == nullptr || outer->containsVertex(bot), "bot not found in graph");
    if(splitGraph == nullptr) {
        throw std::invalid_argument("Cannot call updateGraph until split() has been called");
    }

    outer->removeAllVertices(toSplits.getArraySet());
    std::vector<BaseEdge*> edgesToRemoveVector;
    for(BaseEdge* baseEdge : edgesToRemove) {
        edgesToRemoveVector.emplace_back(baseEdge);
    }
    outer->removeAllEdges(edgesToRemoveVector);

    for(SeqVertex* v : getNewMiddles()) {
        outer->addVertex(v);
    }

    bool hasPrefixSuffixEdge = splitGraph->getEdge(getPrefixV(), getSuffixV()) != nullptr;
    bool hasOnlyPrefixSuffixEdges = hasPrefixSuffixEdge && splitGraph->outDegreeOf(getPrefixV()) == 1;
    bool needPrefixNode = ! getPrefixV()->isEmpty() || (top == nullptr && ! hasOnlyPrefixSuffixEdges);
    bool needSuffixNode = ! getSuffixV()->isEmpty() || (bot == nullptr && ! hasOnlyPrefixSuffixEdges);

    SeqVertex* topForConnect = needPrefixNode ? getPrefixV() : top;
    SeqVertex* botForConnect = needSuffixNode ? getSuffixV() : bot;

    if ( needPrefixNode ) {
        addPrefixNodeAndEdges(top);
    }

    if ( needSuffixNode ) {
        addSuffixNodeAndEdges(bot);
    }

    if ( topForConnect != nullptr ) {
        addEdgesFromTopNode(topForConnect, botForConnect);
    }

    if ( botForConnect != nullptr ) {
        addEdgesToBottomNode(botForConnect);
    }
}

void SharedVertexSequenceSplitter::addPrefixNodeAndEdges(SeqVertex *top) {
    outer->addVertex(getPrefixV());
    if(top != nullptr) {
        outer->addEdge(top, getPrefixV(), BaseEdge::makeOREdge(splitGraph->outgoingEdgesOf(getPrefixV()).getArraySet(), 1));
    }
}

void SharedVertexSequenceSplitter::addSuffixNodeAndEdges(SeqVertex *bot) {
    outer->addVertex(getSuffixV());
    if(bot != nullptr) {
        outer->addEdge(getSuffixV(), bot, BaseEdge::makeOREdge(splitGraph->incomingEdgesOf(getSuffixV()).getArraySet(), 1));
    }
}

void SharedVertexSequenceSplitter::addEdgesFromTopNode(SeqVertex *topForConnect, SeqVertex *botForConnect) {
    for(BaseEdge* e : splitGraph->outgoingEdgesOf(getPrefixV())) {
        SeqVertex* target = splitGraph->getEdgeTarget(e);

        if(target == getSuffixV()) {
            if(botForConnect != nullptr) {
                outer->addEdge(topForConnect, botForConnect, e);
            }
        } else {
            outer->addEdge(topForConnect, target, e);
        }
    }
}

void SharedVertexSequenceSplitter::addEdgesToBottomNode(SeqVertex *botForConnect) {
    for(BaseEdge* e : splitGraph->incomingEdgesOf(getSuffixV())) {
        outer->addEdge(splitGraph->getEdgeSource(e), botForConnect, e);
    }
}
