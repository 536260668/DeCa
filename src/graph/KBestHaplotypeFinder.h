//
// Created by 梦想家xixi on 2021/11/23.
//

#ifndef MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H
#define MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H


#include "SeqGraph.h"
#include "KBestHaplotype.h"
#include <vector>

struct KBestHaplotypeComp {
	bool operator()(const std::shared_ptr<KBestHaplotype> &a, const std::shared_ptr<KBestHaplotype> &b) {
		double aScore = a->getScore(), bScore = b->getScore();
		if (abs(aScore - bScore) > 0.000001)  // a->getScore() != b->getScore();
			return aScore < bScore;
		// convert double to uint64_t
		uint64_t aL, bL;
		memcpy(&aL, &aScore, sizeof(uint64_t));
		memcpy(&bL, &bScore, sizeof(uint64_t));
		if (aL != bL)
			return aL > bL;
		int len1, len2;
		std::shared_ptr<uint8_t[]> base1 = a->getBases(len1);
		std::shared_ptr<uint8_t[]> base2 = b->getBases(len2);
		int len = std::min(len1, len2);
		for (int i = 0; i < len; ++i) {
			if (base1[i] == base2[i]) continue;
			return base1[i] < base2[i];
		}
		return len1 > len2;
	}
};

class KBestHaplotypeFinder {
private:
	std::shared_ptr<SeqGraph> graph;
	std::unordered_set<std::shared_ptr<SeqVertex>> sinks;
	std::unordered_set<std::shared_ptr<SeqVertex>> sources;

	static std::shared_ptr<SeqGraph>
	removeCyclesAndVerticesThatDontLeadToSinks(const std::shared_ptr<SeqGraph> &original,
	                                           std::unordered_set<std::shared_ptr<SeqVertex>> &sources,
	                                           std::unordered_set<std::shared_ptr<SeqVertex>> &sinks);

	static bool findGuiltyVerticesAndEdgesToRemoveCycles(const std::shared_ptr<SeqGraph> &graph,
	                                                     const std::shared_ptr<SeqVertex> &currentVertex,
	                                                     std::unordered_set<std::shared_ptr<SeqVertex>> &sinks,
	                                                     std::unordered_set<std::shared_ptr<BaseEdge>> &edgesToRemove,
	                                                     std::unordered_set<std::shared_ptr<SeqVertex>> &verticesToRemove,
	                                                     std::unordered_set<std::shared_ptr<SeqVertex>> &parentVertices);

public:
	KBestHaplotypeFinder(const std::shared_ptr<SeqGraph> &graph,
	                     std::unordered_set<std::shared_ptr<SeqVertex>> &sources,
	                     std::unordered_set<std::shared_ptr<SeqVertex>> &sinks);

	KBestHaplotypeFinder(std::shared_ptr<SeqGraph> graph, const std::shared_ptr<SeqVertex> &source,
	                     const std::shared_ptr<SeqVertex> &sink);

	KBestHaplotypeFinder(const std::shared_ptr<SeqGraph> &graph);

	std::vector<std::shared_ptr<KBestHaplotype>> findBestHaplotypes(int maxNumberOfHaplotypes);
};


#endif //MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H
