//
// Created by 梦想家xixi on 2021/11/23.
//

#ifndef MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H
#define MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H


#include "SeqGraph.h"
#include "KBestHaplotype.h"
#include <vector>

struct KBestHaplotypeComp{
    bool operator() (std::shared_ptr<KBestHaplotype> a, std::shared_ptr<KBestHaplotype> b) {
        return a->getScore() > b->getScore();
    }
};

class KBestHaplotypeFinder {
private:
    std::shared_ptr<SeqGraph> graph;
    std::unordered_set<std::shared_ptr<SeqVertex>> sinks;
    std::unordered_set<std::shared_ptr<SeqVertex>> sources;
    static std::shared_ptr<SeqGraph> removeCyclesAndVerticesThatDontLeadToSinks(std::shared_ptr<SeqGraph> original, std::unordered_set<std::shared_ptr<SeqVertex>> & sources, std::unordered_set<std::shared_ptr<SeqVertex>> & sinks);
    static bool findGuiltyVerticesAndEdgesToRemoveCycles(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> currentVertex, std::unordered_set<std::shared_ptr<SeqVertex>>& sinks, std::unordered_set<std::shared_ptr<BaseEdge>> & edgesToRemove,
                                                         std::unordered_set<std::shared_ptr<SeqVertex>> & verticesToRemove, std::unordered_set<std::shared_ptr<SeqVertex>> & parentVertices);

public:
    KBestHaplotypeFinder(std::shared_ptr<SeqGraph> graph, std::unordered_set<std::shared_ptr<SeqVertex>> & sources, std::unordered_set<std::shared_ptr<SeqVertex>> & sinks);
    KBestHaplotypeFinder(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex> source, std::shared_ptr<SeqVertex> sink);
    KBestHaplotypeFinder(std::shared_ptr<SeqGraph> graph);
    std::vector<std::shared_ptr<KBestHaplotype>> findBestHaplotypes(int maxNumberOfHaplotypes);
};


#endif //MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H
