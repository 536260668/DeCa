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
private:
    int kmerSize;
public:
    ChainPruner() = default;

    ChainPruner(int kmerSize) : kmerSize(kmerSize) {}

    void pruneLowWeightChains(DirectedSpecifics<V,E> & graph) {
        std::vector<Path<V, E>*> chains = findAllChains(graph);
        std::set<Path<V, E>*> chainsToRemoveset = chainsToRemove(chains);
        for(Path<V, E>* path : chainsToRemoveset)
            graph.removeAllEdges(path->getEdges());
        graph.removeSingletonOrphanVertices();
    }

    std::vector<Path<V, E>*> findAllChains(DirectedSpecifics<V,E> & graph) {
        std::deque<V*> chainStarts;
        std::set<V*> alreadySeen;
        ArraySet<V*> vertexSet = graph.getVertexSet();
        typename std::vector<V*>::iterator viter;
        for(viter = vertexSet.begin(); viter != vertexSet.end(); viter++) {
            if(graph.isSource(*viter)) {
                chainStarts.push_front(*viter);
                alreadySeen.insert(*viter);
            }
        }
        std::vector<Path<V,E>*> chains;
        while(!chainStarts.empty()) {
            V* chainStart = chainStarts.front();
            chainStarts.pop_front();
            ArraySet<E*> outEdges = graph.outgoingEdgesOf(chainStart);
            typename std::vector<E*>::iterator eiter;
            for(eiter = outEdges.begin(); eiter != outEdges.end(); eiter++) {
                Path<V,E>* chain = findChain(*eiter, graph);
                chains.template emplace_back(chain);
                V* chainEnd = chain->getLastVertex();
                if(alreadySeen.find(chainEnd) == alreadySeen.end()) {
                    chainStarts.template emplace_back(chainEnd);
                    alreadySeen.insert(chainEnd);
                }
            }
        }
        return chains;
    }

private:
    Path<V, E>* findChain(E* startEdge, DirectedSpecifics<V,E> & graph) {
        std::vector<E*> edges;
        edges.template emplace_back(startEdge);
        V* firstVertex = graph.getEdgeSource(startEdge);
        V* lastVertex = graph.getEdgeTarget(startEdge);

        while (true) {
            ArraySet<E*> outEdges = graph.outgoingEdgesOf(lastVertex);
            if(outEdges.size() != 1 || graph.inDegreeOf(lastVertex) > 1 || lastVertex == firstVertex) {
                break;
            }
            E* nextEdge = outEdges[0];
            edges.template emplace_back(nextEdge);
            lastVertex = graph.getEdgeTarget(nextEdge);
        }
        return new Path<V,E>(edges, lastVertex, graph, kmerSize);
    }

protected:
    virtual std::set<Path<V,E>*> chainsToRemove(std::vector<Path<V,E>*> chains) = 0;
};


#endif //MUTECT2CPP_MASTER_CHAINPRUNER_H
