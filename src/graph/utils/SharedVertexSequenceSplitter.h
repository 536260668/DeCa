//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_SHAREDVERTEXSEQUENCESPLITTER_H
#define MUTECT2CPP_MASTER_SHAREDVERTEXSEQUENCESPLITTER_H


#include "SeqGraph.h"

class SharedVertexSequenceSplitter {
private:
    SeqGraph* outer;
    SeqVertex* prefixV;
    SeqVertex* suffixV;
    ArraySet<SeqVertex*> toSplits;
    SeqGraph* splitGraph = nullptr;
    std::list<SeqVertex*> newMiddles;
    std::list<BaseEdge*> edgesToRemove;

    BaseEdge* processEdgeToRemove(SeqVertex* v, BaseEdge* e);
    void addPrefixNodeAndEdges(SeqVertex* top);
    void addSuffixNodeAndEdges(SeqVertex* bot);
    void addEdgesFromTopNode(SeqVertex* topForConnect, SeqVertex* botForConnect);
    void addEdgesToBottomNode(SeqVertex* botForConnect);

public:
    SharedVertexSequenceSplitter(SeqGraph* graph, ArraySet<SeqVertex*> toSplitsArg);
    static std::pair<SeqVertex*, SeqVertex*>* commonPrefixAndSuffixOfVertices(ArraySet<SeqVertex*> middleVertices);
    bool meetsMinMergableSequenceForEitherPrefixOrSuffix(int minCommonSequence);
    bool meetsMinMergableSequenceForPrefix(int minCommonSequence);
    bool meetsMinMergableSequenceForSuffix(int minCommonSequence);
    bool splitAndUpdate(SeqVertex* top, SeqVertex* bottom);
    void split();
    void updateGraph(SeqVertex* top, SeqVertex* bot);
    SeqVertex* getPrefixV(){return prefixV;}
    SeqVertex* getSuffixV(){return suffixV;}
    std::list<SeqVertex*> & getNewMiddles() {return newMiddles;}
};


#endif //MUTECT2CPP_MASTER_SHAREDVERTEXSEQUENCESPLITTER_H
