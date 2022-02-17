//
// Created by 梦想家xixi on 2021/11/19.
//

#include "CommonSuffixSplitter.h"
#include <unordered_set>
#include "GraphUtils.h"

bool CommonSuffixSplitter::split(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> v) {
    Mutect2Utils::validateArg(graph.get(), "graph cannot be null");
    Mutect2Utils::validateArg(v.get(), "v cannot be null");
    ArraySet<std::shared_ptr<SeqVertex>> allVertex = graph->getVertexSet();
    Mutect2Utils::validateArg(allVertex.find(v) != allVertex.end(), "graph doesn't contain vertex v ");
    ArraySet<std::shared_ptr<SeqVertex>> toSplit= graph->incomingVerticesOf(v);
    std::shared_ptr<SeqVertex> suffixVTemplate = commonSuffix(graph, v, toSplit);
    if(suffixVTemplate == nullptr) {
        return false;
    }
    std::list<std::shared_ptr<BaseEdge>> edgesToRemove;
    for(std::shared_ptr<SeqVertex> mid : toSplit) {
        std::shared_ptr<SeqVertex> suffixV = std::shared_ptr<SeqVertex>(new SeqVertex(suffixVTemplate->getSequence(), suffixVTemplate->getLength()));
        graph->addVertex(suffixV);
        std::shared_ptr<SeqVertex> prefixV = mid->withoutSuffix(suffixV->getSequence(), suffixV->getLength());
        std::shared_ptr<BaseEdge> out = graph->outgoingEdgeOf(mid);
        std::shared_ptr<SeqVertex> incomingTarget;
        if(prefixV == nullptr) {
            incomingTarget = suffixV;
        } else {
            incomingTarget = prefixV;
            graph->addVertex(prefixV);
            graph->addEdge(prefixV, suffixV, std::shared_ptr<BaseEdge>(new BaseEdge(out->getIsRef(), 1)));
            edgesToRemove.emplace_back(out);
        }
        graph->addEdge(suffixV, graph->getEdgeTarget(out), std::shared_ptr<BaseEdge>(new BaseEdge(out->getIsRef(), out->getMultiplicity())));
        for(std::shared_ptr<BaseEdge> in : graph->incomingEdgesOf(mid)) {
            graph->addEdge(graph->getEdgeSource(in), incomingTarget, std::shared_ptr<BaseEdge>(new BaseEdge(in->getIsRef(), in->getMultiplicity())));
            edgesToRemove.emplace_back(in);
        }
    }
    graph->removeAllVertices(toSplit.getArraySet());
    graph->removeAllEdges(edgesToRemove);
    return true;
}

std::shared_ptr<SeqVertex> CommonSuffixSplitter::commonSuffix(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> v, ArraySet<std::shared_ptr<SeqVertex>> toSplit) {
    if(toSplit.size() < 2) {
        return nullptr;
    } else if (!safeToSplit(graph, v, toSplit)) {
        return nullptr;
    }
    std::shared_ptr<SeqVertex> suffixVTemplate = commonSuffix(toSplit);
    if(suffixVTemplate->isEmpty()) {
        return nullptr;
    } else if ( wouldEliminateRefSource(graph, suffixVTemplate, toSplit) ) {
        return nullptr;
    } else if ( allVerticesAreTheCommonSuffix(suffixVTemplate, toSplit) ) {
        return nullptr;
    } else {
        return suffixVTemplate;
    }
}

bool CommonSuffixSplitter::safeToSplit(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> bot, ArraySet<std::shared_ptr<SeqVertex>> toMerge) {
    ArraySet<std::shared_ptr<SeqVertex>> outgoingVertices = graph->outgoingVerticesOf(bot);
    std::unordered_set<std::shared_ptr<SeqVertex>> outgoingOfBot;
    for(const std::shared_ptr<SeqVertex>& v : outgoingVertices) {
        outgoingOfBot.insert(v);
    }
    for(const std::shared_ptr<SeqVertex>& m : toMerge) {
        ArraySet<std::shared_ptr<BaseEdge>> outs = graph->outgoingEdgesOf(m);
        ArraySet<std::shared_ptr<SeqVertex>> tmp = graph->outgoingVerticesOf(m);
        if ( m == bot || outs.size() != 1 || tmp.find(bot) == tmp.end() ){
            return false;
        }
        if(outgoingOfBot.find(m) != outgoingOfBot.end()) {
            return false;
        }
    }
    return true;
}

std::shared_ptr<SeqVertex> CommonSuffixSplitter::commonSuffix(ArraySet<std::shared_ptr<SeqVertex>> middleVertices) {
    std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> kmers = GraphUtils::getKmers(middleVertices.getArraySet());
    int min = GraphUtils::minKmerLength(kmers);
    int suffixLen = GraphUtils::commonMaximumSuffixLength(kmers, min);
    std::shared_ptr<uint8_t[]> kmer = kmers.begin()->first;
    int kmerLength = kmers.begin()->second;
    int suffixLength;
    std::shared_ptr<uint8_t[]> suffix = Mutect2Utils::copyOfRange(kmer, kmerLength, kmerLength - suffixLen, kmerLength, suffixLength);
    return std::shared_ptr<SeqVertex>(new SeqVertex(suffix, suffixLength));
}

bool
CommonSuffixSplitter::wouldEliminateRefSource(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex>commonSuffix, ArraySet<std::shared_ptr<SeqVertex>> toSplits) {
    for(std::shared_ptr<SeqVertex> toSplit : toSplits) {
        if(graph->isRefSource(toSplit)) {
            return toSplit->getLength() == commonSuffix->getLength();
        }
    }
    return false;
}

bool CommonSuffixSplitter::allVerticesAreTheCommonSuffix(std::shared_ptr<SeqVertex> commonSuffix, ArraySet<std::shared_ptr<SeqVertex>> toSplits) {
    for(std::shared_ptr<SeqVertex> toSplit : toSplits) {
        if(toSplit->getLength() != commonSuffix->getLength()) {
            return false;
        }
    }
    return true;
}
