//
// Created by 梦想家xixi on 2021/11/23.
//

#ifndef MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H
#define MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H


#include "SeqGraph.h"
#include "KBestHaplotype.h"
#include <vector>

struct KBestHaplotypeComp{
    bool operator() (KBestHaplotype* a, KBestHaplotype* b) {
        return a->getScore() > b->getScore();
    }
};

class KBestHaplotypeFinder {
private:
    SeqGraph* graph;
    ArraySet<SeqVertex*> sinks;
    ArraySet<SeqVertex*> sources;
    static SeqGraph* removeCyclesAndVerticesThatDontLeadToSinks(SeqGraph* original, ArraySet<SeqVertex*> & sources, ArraySet<SeqVertex*> & sinks);
    static bool findGuiltyVerticesAndEdgesToRemoveCycles(SeqGraph* graph, SeqVertex* currentVertex, ArraySet<SeqVertex*>& sinks, std::set<BaseEdge*> & edgesToRemove, std::set<SeqVertex*> & verticesToRemove, std::set<SeqVertex*> & parentVertices);

public:
    KBestHaplotypeFinder(SeqGraph* graph, ArraySet<SeqVertex*> & sources, ArraySet<SeqVertex*> & sinks);
    KBestHaplotypeFinder(SeqGraph* graph, SeqVertex* source, SeqVertex* sink);
    KBestHaplotypeFinder(SeqGraph* graph);
    std::vector<KBestHaplotype*> findBestHaplotypes(int maxNumberOfHaplotypes);
};


#endif //MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H
