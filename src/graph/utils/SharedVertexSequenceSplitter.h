//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_SHAREDVERTEXSEQUENCESPLITTER_H
#define MUTECT2CPP_MASTER_SHAREDVERTEXSEQUENCESPLITTER_H


#include "SeqGraph.h"

class SharedVertexSequenceSplitter {
private:
    std::shared_ptr<SeqGraph> outer;
    std::shared_ptr<SeqVertex> prefixV;
    std::shared_ptr<SeqVertex> suffixV;
    std::unordered_set<std::shared_ptr<SeqVertex>> toSplits;
    std::shared_ptr<SeqGraph> splitGraph = nullptr;
    std::list<std::shared_ptr<SeqVertex>> newMiddles;
    std::list<std::shared_ptr<BaseEdge>> edgesToRemove;

    std::shared_ptr<BaseEdge> processEdgeToRemove(std::shared_ptr<SeqVertex> v, std::shared_ptr<BaseEdge> e);
    void addPrefixNodeAndEdges(std::shared_ptr<SeqVertex> top);
    void addSuffixNodeAndEdges(std::shared_ptr<SeqVertex> bot);
    void addEdgesFromTopNode(std::shared_ptr<SeqVertex> topForConnect, std::shared_ptr<SeqVertex> botForConnect);
    void addEdgesToBottomNode(std::shared_ptr<SeqVertex> botForConnect);

public:
    SharedVertexSequenceSplitter(std::shared_ptr<SeqGraph> graph, std::unordered_set<std::shared_ptr<SeqVertex>> toSplitsArg);
    static std::pair<std::shared_ptr<SeqVertex>, std::shared_ptr<SeqVertex>> commonPrefixAndSuffixOfVertices(std::unordered_set<std::shared_ptr<SeqVertex>> middleVertices);
    bool meetsMinMergableSequenceForEitherPrefixOrSuffix(int minCommonSequence);
    bool meetsMinMergableSequenceForPrefix(int minCommonSequence);
    bool meetsMinMergableSequenceForSuffix(int minCommonSequence);
    bool splitAndUpdate(std::shared_ptr<SeqVertex> top, std::shared_ptr<SeqVertex> bottom);
    void split();
    void updateGraph(std::shared_ptr<SeqVertex> top, std::shared_ptr<SeqVertex> bot);
    std::shared_ptr<SeqVertex> getPrefixV(){return prefixV;}
    std::shared_ptr<SeqVertex> getSuffixV(){return suffixV;}
    std::list<std::shared_ptr<SeqVertex>> & getNewMiddles() {return newMiddles;}
};


#endif //MUTECT2CPP_MASTER_SHAREDVERTEXSEQUENCESPLITTER_H
