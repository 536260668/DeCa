//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_GRAPHUTILS_H
#define MUTECT2CPP_MASTER_GRAPHUTILS_H

#include <list>
#include <set>
#include "SeqVertex.h"
#include "BaseGraph/DirectedSpecifics.h"

class GraphUtils {
public:
    static int minKmerLength(std::list<std::pair<uint8_t *, int>> & kmers);

    static int commonMaximumPrefixLength(std::list<std::pair<uint8_t *, int>> & kmers);

    static int commonMaximumSuffixLength(std::list<std::pair<uint8_t *, int>> & listOfBytes, int minLength);

    static std::list<std::pair<uint8_t *, int>> getKmers(std::vector<SeqVertex*> vertices);

    template<class T, class E>
            static bool graphEquals(DirectedSpecifics<T,E>* g1, DirectedSpecifics<T,E>* g2){
        Mutect2Utils::validateArg(g1, "g1");
        Mutect2Utils::validateArg(g2, "g2");
        ArraySet<T*> vertices1 = g1->getVertexSet();
        ArraySet<T*> vertices2 = g2->getVertexSet();
        ArraySet<E*> edges1 = g1->getEdgeSet();
        ArraySet<E*> edges2 = g2->getEdgeSet();
        if(vertices1.size() != vertices2.size() || edges1.size() != edges2.size())
            return false;
        return true;
    };
};


#endif //MUTECT2CPP_MASTER_GRAPHUTILS_H
