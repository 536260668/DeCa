//
// Created by 梦想家xixi on 2021/10/27.
//

#ifndef MUTECT2CPP_MASTER_PATH_H
#define MUTECT2CPP_MASTER_PATH_H

#include <cstring>
#include <vector>
#include "graph/BaseGraph/DirectedSpecifics.h"
#include "Mutect2Utils.h"
#include <iostream>

template<class T, class E>
class Path {
private:
    // the last vertex seen in the path
    std::shared_ptr<T> lastVertex;

    // the list of edges comprising the path
    std::vector<std::shared_ptr<E>> edgesInOrder;

    // the graph from which this path originated
    std::shared_ptr<DirectedSpecifics<T, E>>  graph;


public:
    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    Path(std::shared_ptr<T> initialVertex, std::shared_ptr<DirectedSpecifics<T, E>> graph) :  lastVertex(initialVertex), graph(graph){
        Mutect2Utils::validateArg(initialVertex.get(), "initialVertex cannot be null");
        Mutect2Utils::validateArg(graph->containsVertex(initialVertex), "Vertex must be part of graph.");
    }

    /**
     * Constructor that does not check arguments' validity i.e. doesn't check that edges are in order
     */
    Path(std::vector<std::shared_ptr<E>> edgesInOrder, std::shared_ptr<T> lastVertex, std::shared_ptr<DirectedSpecifics<T, E>> graph) : lastVertex(lastVertex), graph(graph), edgesInOrder(edgesInOrder){}

    /**
     * Create a new Path extending p with edge
     *
     * @param p the path to extend.
     * @param edge the edge to extend path with.
     *
     * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code p}'s graph, or {@code edge} does not have as a source the last vertex in {@code p}.
     */
    Path(Path<T, E> p, std::shared_ptr<E> edge) : graph(p.graph), lastVertex(p.graph->getEdgeTarget(edge)) {
        Mutect2Utils::validateArg(edge.get(), "Edge cannot be null");
        Mutect2Utils::validateArg(p.graph->getEdgeSource(edge) == p.lastVertex, "Edges added to path must be contiguous.");
        for(typename std::vector<std::shared_ptr<E>>::iterator iter = p.edgesInOrder.begin(); iter != p.edgesInOrder.end(); iter++) {
            edgesInOrder.template emplace_back(*iter);
        }
        edgesInOrder.template emplace_back(edge);
    }

    int length() const {return edgesInOrder.size();}

    Path(std::shared_ptr<E> edge, Path<T,E> p) : graph(p.graph), lastVertex(p.lastVertex){
         Mutect2Utils::validateArg(edge, "Edge cannot be null");
         Mutect2Utils::validateArg(p.graph.containsEdge(edge), "Graph must contain edge ");
         Mutect2Utils::validateArg(p.graph.getEdgeTarget(edge) == p.getFirstVertex(), "Edges added to path must be contiguous");
        for(typename std::vector<std::shared_ptr<E>>::iterator iter = p.edgesInOrder.begin(); iter != p.edgesInOrder.end(); iter++) {
            edgesInOrder.insert(*iter);
        }
        edgesInOrder.insert(edge);
     }

     bool containsVertex(std::shared_ptr<T> v) {
         Mutect2Utils::validateArg(v, "Vertex cannot be null");
         std::vector<std::shared_ptr<T>> res;
         res.emplace_back(getFirstVertex());
         for(typename std::vector<std::shared_ptr<E>>::iterator iter = edgesInOrder.begin(); iter != edgesInOrder.end(); iter++) {
             res.emplace_back(graph->getEdgeTarget(*iter));
         }
         return std::find(res.begin(), res.end(), v) != res.end();
     }

    std::shared_ptr<DirectedSpecifics<T, E>> getGraph() {return graph;}

    std::vector<std::shared_ptr<E>> & getEdges() {return edgesInOrder;}

    std::shared_ptr<E> getLastEdge() {return edgesInOrder[edgesInOrder.size()-1];}

    std::vector<std::shared_ptr<T>> getVertices() {
         std::vector<std::shared_ptr<T>> res;
         res.emplace_back(getFirstVertex());
         for(typename std::vector<std::shared_ptr<E>>::iterator iter = edgesInOrder.begin(); iter != edgesInOrder.end(); iter++) {
             res.emplace_back(graph->getEdgeTarget(*iter));
         }
         return res;
     }

    std::shared_ptr<T> getFirstVertex() {
         if(edgesInOrder.empty()) {
             return lastVertex;
         } else {
             return getGraph()->getEdgeSource(edgesInOrder[0]);
         }
     }

     //返回指针需要delete
     std::shared_ptr<uint8_t[]> getBases(int & reslength) {
         if(getEdges().empty()) {
             return graph->getAdditionalSequence(lastVertex);
         }
         std::shared_ptr<T> source = graph->getEdgeSource(edgesInOrder[0]);
         std::shared_ptr<uint8_t[]> bases = graph->getAdditionalSequence(source);
         int basesLength = graph->getAdditionalSequenceLength(source);
         std::shared_ptr<uint8_t[]> res(new uint8_t[600]);
         memcpy(res.get(), bases.get(), basesLength);
         int length = basesLength;
         int start = basesLength;
         for(int i = 0; i < edgesInOrder.size(); i++) {
             std::shared_ptr<T> target = graph->getEdgeTarget(edgesInOrder[i]);
//             if(length <= start) {
//                 length *= 2;
//                 std::shared_ptr<uint8_t[]> tmp(new uint8_t[length]);
//                 memcpy(tmp.get(), res.get(), start);
//                 res = tmp;
//             }
             bases = graph->getAdditionalSequence(target);
             basesLength = graph->getAdditionalSequenceLength(target);
             memcpy(res.get()+start, bases.get(), basesLength);
             start += basesLength;
         }
         std::shared_ptr<uint8_t[]> tmp1(new uint8_t[start]);
         memcpy(tmp1.get(), res.get(), start);
//         std::shared_ptr<uint8_t> test = res;
         reslength = start;
//         std::cout << start << std::endl;
//         for(int i = 0; i < reslength; i++) {
//             std::cout << tmp.get()[i];
//         }
//         std::cout << std::endl;
         return tmp1;
     }

    std::shared_ptr<T> getLastVertex() {return lastVertex;}
};

#endif //MUTECT2CPP_MASTER_PATH_H
