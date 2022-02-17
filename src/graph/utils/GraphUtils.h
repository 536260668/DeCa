//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_GRAPHUTILS_H
#define MUTECT2CPP_MASTER_GRAPHUTILS_H

#include <list>
#include <set>
#include "SeqVertex.h"
#include "BaseEdge.h"
#include "BaseGraph/DirectedSpecifics.h"

class GraphUtils {
public:
    static int minKmerLength(std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> & kmers);

    static int commonMaximumPrefixLength(std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> & kmers);

    static int commonMaximumSuffixLength(std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> & listOfBytes, int minLength);

    static std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> getKmers(std::vector<std::shared_ptr<SeqVertex>> vertices);

    template<class T, class E>
            static bool graphEquals(DirectedSpecifics<T,E>* g1, DirectedSpecifics<T,E>* g2){
        Mutect2Utils::validateArg(g1, "g1");
        Mutect2Utils::validateArg(g2, "g2");
        ArraySet<std::shared_ptr<T>> vertices1 = g1->getVertexSet();
        ArraySet<std::shared_ptr<T>> vertices2 = g2->getVertexSet();
        ArraySet<std::shared_ptr<E>> edges1 = g1->getEdgeSet();
        ArraySet<std::shared_ptr<E>> edges2 = g2->getEdgeSet();
        if(vertices1.size() != vertices2.size() || edges1.size() != edges2.size())
            return false;
        for(std::shared_ptr<BaseVertex> v1 : vertices1) {
            bool flag = false;
            for(std::shared_ptr<BaseVertex> v2 : vertices2) {
                if(*v1 == *v2) {
                    flag = true;
                    break;
                }
            }
            if(!flag) {
                return false;
            }
        }
        for(std::shared_ptr<BaseEdge> edge1 : edges1) {
            bool flag = false;
            for(std::shared_ptr<BaseEdge> edge2 : edges2) {
                if(*g1->getEdgeSource(edge1) == *g2->getEdgeSource(edge2) && *g1->getEdgeTarget(edge1) == *g2->getEdgeTarget(edge2)) {
                    flag = true;
                    break;
                }
            }
            if(!flag) {
                return false;
            }
        }
        return true;
    };
};


#endif //MUTECT2CPP_MASTER_GRAPHUTILS_H
