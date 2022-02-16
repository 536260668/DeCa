//
// Created by 梦想家xixi on 2021/10/28.
//

#ifndef MUTECT2CPP_MASTER_CHAINPRUNER_H
#define MUTECT2CPP_MASTER_CHAINPRUNER_H

#include <vector>
#include "./Path.h"
#include <deque>
#include <set>
#include "graph/set/ArraySet.h"


template<class V, class E>
class ChainPruner {
public:
    ChainPruner() = default;


    void pruneLowWeightChains(DirectedSpecifics<V,E> & graph) {
        std::vector<Path<V, E>*> chains = findAllChains(graph);
        std::set<Path<V, E>*> chainsToRemoveset = chainsToRemove(chains);
        for(Path<V, E>* path : chainsToRemoveset)
            graph.removeAllEdges(path->getEdges());
        graph.removeSingletonOrphanVertices();
    }

    std::vector<Path<V, E>*> findAllChains(DirectedSpecifics<V,E> & graph) {
        std::deque<std::shared_ptr<V>> chainStarts;
        std::set<std::shared_ptr<V>> alreadySeen;
        ArraySet<std::shared_ptr<V>> vertexSet = graph.getVertexSet();
        typename std::vector<std::shared_ptr<V>>::iterator viter;
        for(viter = vertexSet.begin(); viter != vertexSet.end(); viter++) {
            if(graph.isSource(*viter)) {
                chainStarts.push_front(*viter);
                alreadySeen.insert(*viter);
            }
        }
        std::vector<Path<V,E>*> chains;
        while(!chainStarts.empty()) {
            std::shared_ptr<V> chainStart = chainStarts.front();
            chainStarts.pop_front();
            ArraySet<std::shared_ptr<E>> outEdges = graph.outgoingEdgesOf(chainStart);
            typename std::vector<std::shared_ptr<E>>::iterator eiter;
            for(eiter = outEdges.begin(); eiter != outEdges.end(); eiter++) {
                Path<V,E>* chain = findChain(*eiter, graph);
                chains.template emplace_back(chain);
                std::shared_ptr<V> chainEnd = chain->getLastVertex();
                if(alreadySeen.find(chainEnd) == alreadySeen.end()) {
                    chainStarts.template emplace_back(chainEnd);
                    alreadySeen.insert(chainEnd);
                }
            }
        }
        return chains;
    }

private:
    Path<V, E>* findChain(std::shared_ptr<E> startEdge, DirectedSpecifics<V,E> & graph) {
        std::vector<std::shared_ptr<E>> edges;
        edges.template emplace_back(startEdge);
        std::shared_ptr<V> firstVertex = graph.getEdgeSource(startEdge);
        std::shared_ptr<V> lastVertex = graph.getEdgeTarget(startEdge);

        while (true) {
            ArraySet<std::shared_ptr<E>> outEdges = graph.outgoingEdgesOf(lastVertex);
            if(outEdges.size() != 1 || graph.inDegreeOf(lastVertex) > 1 || lastVertex == firstVertex) {
                break;
            }
            std::shared_ptr<E> nextEdge = outEdges[0];
            edges.template emplace_back(nextEdge);
            lastVertex = graph.getEdgeTarget(nextEdge);
        }
        return new Path<V,E>(edges, lastVertex, graph);
    }

protected:
    virtual std::set<Path<V,E>*> chainsToRemove(std::vector<Path<V,E>*> chains) = 0;
};


#endif //MUTECT2CPP_MASTER_CHAINPRUNER_H
