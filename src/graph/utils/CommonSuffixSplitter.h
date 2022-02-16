//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_COMMONSUFFIXSPLITTER_H
#define MUTECT2CPP_MASTER_COMMONSUFFIXSPLITTER_H


#include "SeqGraph.h"

class CommonSuffixSplitter {
public:
    static bool split(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> v);

private:
    static std::shared_ptr<SeqVertex> commonSuffix(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> v, ArraySet<std::shared_ptr<SeqVertex>> toSplit);
    static bool safeToSplit(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> bot, ArraySet<std::shared_ptr<SeqVertex>> toSplit);
    static std::shared_ptr<SeqVertex> commonSuffix(ArraySet<std::shared_ptr<SeqVertex>> toSplit);
    static bool wouldEliminateRefSource(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> commonSuffix, ArraySet<std::shared_ptr<SeqVertex>> toSplit);
    static bool allVerticesAreTheCommonSuffix(std::shared_ptr<SeqVertex> commonSuffix, ArraySet<std::shared_ptr<SeqVertex>> toSplits);
};


#endif //MUTECT2CPP_MASTER_COMMONSUFFIXSPLITTER_H
