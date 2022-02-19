//
// Created by 梦想家xixi on 2021/11/19.
//

#include "SharedVertexSequenceSplitter.h"
#include "GraphUtils.h"

std::pair<std::shared_ptr<SeqVertex> , std::shared_ptr<SeqVertex>>
SharedVertexSequenceSplitter::commonPrefixAndSuffixOfVertices(std::set<std::shared_ptr<SeqVertex>> middleVertices){
    std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> kmers;
    int min = INT32_MAX;
    for(std::shared_ptr<SeqVertex> v : middleVertices) {
        std::pair<std::shared_ptr<uint8_t[]>, int> tmp;
        tmp.first = v->getSequence();
        tmp.second = v->getLength();
        kmers.emplace_back(tmp);
        min = std::min(min, tmp.second);
    }

    int prefixLen = GraphUtils::commonMaximumPrefixLength(kmers);
    int suffixLen = GraphUtils::commonMaximumSuffixLength(kmers, min - prefixLen);

    std::shared_ptr<uint8_t[]> kmer = kmers.begin()->first;
    int length = kmers.begin()->second;
    int prefixLength;
    std::shared_ptr<uint8_t[]> prefix = Mutect2Utils::copyOfRange(kmer, length, 0, prefixLen, prefixLength);
    int suffixLength;
    std::shared_ptr<uint8_t[]> suffix = Mutect2Utils::copyOfRange(kmer, length, length - suffixLen, length, suffixLength);
    return std::pair<std::shared_ptr<SeqVertex>, std::shared_ptr<SeqVertex>>(new SeqVertex(prefix, prefixLength), new SeqVertex(suffix, suffixLength));
}

SharedVertexSequenceSplitter::SharedVertexSequenceSplitter(std::shared_ptr<SeqGraph> graph, std::set<std::shared_ptr<SeqVertex>> toSplitsArg) : outer(graph), toSplits(toSplitsArg){
    Mutect2Utils::validateArg(graph.get(), "graph cannot be null");
    Mutect2Utils::validateArg(toSplitsArg.size() > 1, "Can only split at least 2 vertices");
    for(std::shared_ptr<SeqVertex> v : toSplitsArg) {
        std::set<std::shared_ptr<SeqVertex>>& allVertex = graph->getVertexSet();
        if(allVertex.find(v) == allVertex.end())
            throw std::invalid_argument("graph doesn't contain all of the vertices to split");
    }
    std::pair<std::shared_ptr<SeqVertex>, std::shared_ptr<SeqVertex>> prefixAndSuffix = commonPrefixAndSuffixOfVertices(toSplits);
    prefixV = prefixAndSuffix.first;
    suffixV = prefixAndSuffix.second;
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

bool SharedVertexSequenceSplitter::splitAndUpdate(std::shared_ptr<SeqVertex> top, std::shared_ptr<SeqVertex> bottom) {
    split();
    updateGraph(top, bottom);
    return true;
}

void SharedVertexSequenceSplitter::split() {
    splitGraph = std::shared_ptr<SeqGraph>(new SeqGraph(outer->getKmerSize()));
    splitGraph->addVertex(getPrefixV());
    splitGraph->addVertex(getSuffixV());
    for(std::shared_ptr<SeqVertex> mid : toSplits) {
        std::shared_ptr<BaseEdge> toMid = processEdgeToRemove(mid, outer->incomingEdgeOf(mid));
        std::shared_ptr<BaseEdge> fromMid = processEdgeToRemove(mid, outer->outgoingEdgeOf(mid));

        std::shared_ptr<SeqVertex> remaining = mid->withoutPrefixAndSuffix(getPrefixV()->getSequence(), getPrefixV()->getLength(), getSuffixV()->getSequence(), getSuffixV()->getLength());
        if(remaining != nullptr) {
            splitGraph->addVertex(remaining);
            getNewMiddles().emplace_back(remaining);
            splitGraph->addEdge(getPrefixV(), remaining, toMid);
            splitGraph->addEdge(remaining, getSuffixV(), fromMid);
        } else {
            std::shared_ptr<BaseEdge> tmp(new BaseEdge(toMid->getIsRef(), toMid->getMultiplicity()));
            tmp->add(*fromMid);
            splitGraph->addOrUpdateEdge(getPrefixV(), getSuffixV(), tmp);
        }
    }
}

std::shared_ptr<BaseEdge> SharedVertexSequenceSplitter::processEdgeToRemove(std::shared_ptr<SeqVertex> v, std::shared_ptr<BaseEdge> e) {
    if (e == nullptr) {
        return std::shared_ptr<BaseEdge>(new BaseEdge(outer->isReferenceNode(v), 0));
    } else {
        edgesToRemove.emplace_back(e);
        return std::shared_ptr<BaseEdge>(new BaseEdge(e->getIsRef(), e->getMultiplicity()));
    }
}

void SharedVertexSequenceSplitter::updateGraph(std::shared_ptr<SeqVertex> top, std::shared_ptr<SeqVertex> bot) {
    for(std::shared_ptr<SeqVertex> v : toSplits) {
        std::set<std::shared_ptr<SeqVertex>>& allVertex = outer->getVertexSet();
        if(allVertex.find(v) == allVertex.end())
            throw std::invalid_argument("graph doesn't contain all of the vertices to split");
    }
    Mutect2Utils::validateArg(top != nullptr || bot != nullptr, "Cannot update graph without at least one top or bot vertex, but both were null");
    Mutect2Utils::validateArg(top == nullptr || outer->containsVertex(top), "top not found in graph");
    Mutect2Utils::validateArg(bot == nullptr || outer->containsVertex(bot), "bot not found in graph");
    if(splitGraph == nullptr) {
        throw std::invalid_argument("Cannot call updateGraph until split() has been called");
    }

    outer->removeAllVertices(toSplits);
    std::vector<std::shared_ptr<BaseEdge>> edgesToRemoveVector;
    for(std::shared_ptr<BaseEdge> baseEdge : edgesToRemove) {
        edgesToRemoveVector.emplace_back(baseEdge);
    }
    outer->removeAllEdges(edgesToRemoveVector);

    for(std::shared_ptr<SeqVertex> v : getNewMiddles()) {
        outer->addVertex(v);
    }

    bool hasPrefixSuffixEdge = splitGraph->getEdge(getPrefixV(), getSuffixV()) != nullptr;
    bool hasOnlyPrefixSuffixEdges = hasPrefixSuffixEdge && splitGraph->outDegreeOf(getPrefixV()) == 1;
    bool needPrefixNode = ! getPrefixV()->isEmpty() || (top == nullptr && ! hasOnlyPrefixSuffixEdges);
    bool needSuffixNode = ! getSuffixV()->isEmpty() || (bot == nullptr && ! hasOnlyPrefixSuffixEdges);

    std::shared_ptr<SeqVertex> topForConnect = needPrefixNode ? getPrefixV() : top;
    std::shared_ptr<SeqVertex> botForConnect = needSuffixNode ? getSuffixV() : bot;

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

void SharedVertexSequenceSplitter::addPrefixNodeAndEdges(std::shared_ptr<SeqVertex> top) {
    outer->addVertex(getPrefixV());
    if(top != nullptr) {
        outer->addEdge(top, getPrefixV(), BaseEdge::makeOREdge(splitGraph->outgoingEdgesOf(getPrefixV()), 1));
    }
}

void SharedVertexSequenceSplitter::addSuffixNodeAndEdges(std::shared_ptr<SeqVertex> bot) {
    outer->addVertex(getSuffixV());
    if(bot != nullptr) {
        outer->addEdge(getSuffixV(), bot, BaseEdge::makeOREdge(splitGraph->incomingEdgesOf(getSuffixV()), 1));
    }
}

void SharedVertexSequenceSplitter::addEdgesFromTopNode(std::shared_ptr<SeqVertex> topForConnect, std::shared_ptr<SeqVertex> botForConnect) {
    for(std::shared_ptr<BaseEdge> e : splitGraph->outgoingEdgesOf(getPrefixV())) {
        std::shared_ptr<SeqVertex> target = splitGraph->getEdgeTarget(e);

        if(target == getSuffixV()) {
            if(botForConnect != nullptr) {
                outer->addEdge(topForConnect, botForConnect, e);
            }
        } else {
            outer->addEdge(topForConnect, target, e);
        }
    }
}

void SharedVertexSequenceSplitter::addEdgesToBottomNode(std::shared_ptr<SeqVertex> botForConnect) {
    for(std::shared_ptr<BaseEdge> e : splitGraph->incomingEdgesOf(getSuffixV())) {
        outer->addEdge(splitGraph->getEdgeSource(e), botForConnect, e);
    }
}
