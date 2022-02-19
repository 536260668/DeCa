//
// Created by 梦想家xixi on 2021/11/20.
//

#ifndef MUTECT2CPP_MASTER_SHAREDSEQUENCEMERGER_H
#define MUTECT2CPP_MASTER_SHAREDSEQUENCEMERGER_H


#include "SeqGraph.h"

class SharedSequenceMerger {
public:
    static bool canMerge(const std::shared_ptr<SeqGraph>& graph, std::shared_ptr<SeqVertex> v, std::set<std::shared_ptr<SeqVertex>> incomingVertices);
    static bool merge(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> v);
};


#endif //MUTECT2CPP_MASTER_SHAREDSEQUENCEMERGER_H
