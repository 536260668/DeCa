//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_COMMONSUFFIXSPLITTER_H
#define MUTECT2CPP_MASTER_COMMONSUFFIXSPLITTER_H


#include "SeqGraph.h"

class CommonSuffixSplitter {
public:
    static bool split(SeqGraph* graph, SeqVertex* v);

private:
    static SeqVertex* commonSuffix(SeqGraph* graph, SeqVertex* v, ArraySet<SeqVertex*> toSplit);
    static bool safeToSplit(SeqGraph* graph, SeqVertex* bot, ArraySet<SeqVertex*> toSplit);
    static SeqVertex* commonSuffix(ArraySet<SeqVertex*> toSplit);
    static bool wouldEliminateRefSource(SeqGraph* graph, SeqVertex* commonSuffix, ArraySet<SeqVertex*> toSplit);
    static bool allVerticesAreTheCommonSuffix(SeqVertex* commonSuffix, ArraySet<SeqVertex*> toSplits);
};


#endif //MUTECT2CPP_MASTER_COMMONSUFFIXSPLITTER_H
