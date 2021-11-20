//
// Created by 梦想家xixi on 2021/11/19.
//

#include "CommonSuffixSplitter.h"
#include <unordered_set>
#include "GraphUtils.h"

bool CommonSuffixSplitter::split(SeqGraph *graph, SeqVertex *v) {
    Mutect2Utils::validateArg(graph, "graph cannot be null");
    Mutect2Utils::validateArg(v, "v cannot be null");
    ArraySet<SeqVertex*> allVertex = graph->getVertexSet();
    Mutect2Utils::validateArg(allVertex.find(v) != allVertex.end(), "graph doesn't contain vertex v ");
    ArraySet<SeqVertex*> toSplit= graph->incomingVerticesOf(v);
    SeqVertex* suffixVTemplate = commonSuffix(graph, v, toSplit);
    if(suffixVTemplate == nullptr) {
        return false;
    }
    std::list<BaseEdge*> edgesToRemove;
    for(SeqVertex* mid : toSplit) {
        SeqVertex* suffixV = new SeqVertex(suffixVTemplate->getSequence(), suffixVTemplate->getLength());
        graph->addVertex(suffixV);
        SeqVertex* prefixV = mid->withoutSuffix(suffixV->getSequence(), suffixV->getLength());
        BaseEdge* out = graph->outgoingEdgeOf(mid);
        SeqVertex* incomingTarget;
        if(prefixV == nullptr) {
            incomingTarget = suffixV;
        } else {
            incomingTarget = prefixV;
            graph->addVertex(prefixV);
            graph->addEdge(prefixV, suffixV, new BaseEdge(out->getIsRef(), 1));
            edgesToRemove.emplace_back(out);
        }
        graph->addEdge(suffixV, graph->getEdgeTarget(out), new BaseEdge(out->getIsRef(), out->getMultiplicity()));
        for(BaseEdge* in : graph->incomingEdgesOf(mid)) {
            graph->addEdge(graph->getEdgeSource(in), incomingTarget, new BaseEdge(in->getIsRef(), in->getMultiplicity()));
            edgesToRemove.emplace_back(in);
        }
    }
    graph->removeAllVertices(toSplit.getArraySet());
    graph->removeAllEdges(std::vector<BaseEdge*>(edgesToRemove.begin(), edgesToRemove.end()));
    return true;
}

SeqVertex* CommonSuffixSplitter::commonSuffix(SeqGraph *graph, SeqVertex *v, ArraySet<SeqVertex *> toSplit) {
    if(toSplit.size() < 2) {
        return nullptr;
    } else if (!safeToSplit(graph, v, toSplit)) {
        return nullptr;
    }
    SeqVertex* suffixVTemplate = commonSuffix(toSplit);
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

bool CommonSuffixSplitter::safeToSplit(SeqGraph *graph, SeqVertex* bot, ArraySet<SeqVertex *> toMerge) {
    ArraySet<SeqVertex*> outgoingVertices = graph->outgoingVerticesOf(bot);
    std::unordered_set<SeqVertex*> outgoingOfBot;
    for(SeqVertex* v : outgoingVertices) {
        outgoingOfBot.insert(v);
    }
    for(SeqVertex* m : toMerge) {
        ArraySet<BaseEdge*> outs = graph->outgoingEdgesOf(m);
        ArraySet<SeqVertex*> tmp = graph->outgoingVerticesOf(m);
        if ( m == bot || outs.size() != 1 || tmp.find(bot) == tmp.end() ){
            return false;
        }
        if(outgoingOfBot.find(m) != outgoingOfBot.end()) {
            return false;
        }
    }
    return true;
}

SeqVertex *CommonSuffixSplitter::commonSuffix(ArraySet<SeqVertex *> middleVertices) {
    std::list<std::pair<uint8_t *, int>> kmers = GraphUtils::getKmers(middleVertices.getArraySet());
    int min = GraphUtils::minKmerLength(kmers);
    int suffixLen = GraphUtils::commonMaximumSuffixLength(kmers, min);
    uint8_t * kmer = kmers.begin()->first;
    int kmerLength = kmers.begin()->second;
    int suffixLength;
    uint8_t * suffix = Mutect2Utils::copyOfRange(kmer, kmerLength, kmerLength - suffixLen, kmerLength, suffixLength);
    return new SeqVertex(suffix, suffixLength);
}

bool
CommonSuffixSplitter::wouldEliminateRefSource(SeqGraph *graph, SeqVertex *commonSuffix, ArraySet<SeqVertex *> toSplits) {
    for(SeqVertex* toSplit : toSplits) {
        if(graph->isRefSource(toSplit)) {
            return toSplit->getLength() == commonSuffix->getLength();
        }
    }
    return false;
}

bool CommonSuffixSplitter::allVerticesAreTheCommonSuffix(SeqVertex *commonSuffix, ArraySet<SeqVertex *> toSplits) {
    for(SeqVertex* toSplit : toSplits) {
        if(toSplit->getLength() != commonSuffix->getLength()) {
            return false;
        }
    }
    return true;
}
