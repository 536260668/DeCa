//
// Created by 梦想家xixi on 2021/10/27.
//

#ifndef MUTECT2CPP_MASTER_PATH_H
#define MUTECT2CPP_MASTER_PATH_H

#include <cstring>
#include <vector>
#include "graph/BaseGraph/DirectedSpecifics.h"
#include "Mutect2Utils.h"

template<class T, class E>
class Path {
private:
    // the last vertex seen in the path
    T* lastVertex;

    // the list of edges comprising the path
    std::vector<E*> edgesInOrder;

    // the graph from which this path originated
    DirectedSpecifics<T, E> graph;

    int kmerSize;

public:
    /**
     * Create a new Path containing no edges and starting at initialVertex
     * @param initialVertex the starting vertex of the path
     * @param graph the graph this path will follow through
     */
    Path(T* initialVertex, DirectedSpecifics<T, E> graph, int kmerSize) :  lastVertex(lastVertex), graph(graph), kmerSize(kmerSize){
        Mutect2Utils::validateArg(initialVertex, "initialVertex cannot be null");
        Mutect2Utils::validateArg(graph.containsVertex(initialVertex), "Vertex must be part of graph.");
    }

    /**
     * Constructor that does not check arguments' validity i.e. doesn't check that edges are in order
     */
    Path(std::vector<E*> edgesInOrder, T* lastVertex, DirectedSpecifics<T, E> graph, int kmerSize) : lastVertex(lastVertex), graph(graph), edgesInOrder(edgesInOrder), kmerSize(kmerSize){}

    /**
     * Create a new Path extending p with edge
     *
     * @param p the path to extend.
     * @param edge the edge to extend path with.
     *
     * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
     * not part of {@code p}'s graph, or {@code edge} does not have as a source the last vertex in {@code p}.
     */
     Path(Path<T, E> p, E* edge, int kmerSize) : graph(p.graph), lastVertex(p.graph.getEdgeTarget(edge)), kmerSize(kmerSize) {
         Mutect2Utils::validateArg(edge, "Edge cannot be null");
         Mutect2Utils::validateArg(p.graph.getEdgeTarget(edge) == p.lastVertex, "Edges added to path must be contiguous.");
         for(typename std::vector<E*>::iterator iter = p.edgesInOrder.begin(); iter != p.edgesInOrder.end(); iter++) {
             edgesInOrder.insert(*iter);
         }
         edgesInOrder.insert(edge);
     }

    int length() const {return edgesInOrder.size();}

    Path(E* edge, Path<T,E> p, int kmerSize) : graph(p.graph), lastVertex(p.lastVertex), kmerSize(kmerSize){
         Mutect2Utils::validateArg(edge, "Edge cannot be null");
         Mutect2Utils::validateArg(p.graph.containsEdge(edge), "Graph must contain edge ");
         Mutect2Utils::validateArg(p.graph.getEdgeTarget(edge) == p.getFirstVertex(), "Edges added to path must be contiguous");
        for(typename std::vector<E*>::iterator iter = p.edgesInOrder.begin(); iter != p.edgesInOrder.end(); iter++) {
            edgesInOrder.insert(*iter);
        }
        edgesInOrder.insert(edge);
     }

     bool containsVertex(T* v) {
         Mutect2Utils::validateArg(v, "Vertex cannot be null");
         std::vector<T*> res;
         res.emplace_back(getFirstVertex());
         for(typename std::vector<E*>::iterator iter = edgesInOrder.begin(); iter != edgesInOrder.end(); iter++) {
             res.emplace_back(graph.getEdgeTarget(*iter));
         }
         return std::find(res.begin(), res.end(), v) != res.end();
     }

    DirectedSpecifics<T, E> getGraph() {return graph;}

    std::vector<E*> & getEdges() {return edgesInOrder;}

    E* getLastEdge() {return edgesInOrder[edgesInOrder.size()-1];}

    T* getLastVertex() {return lastVertex;}

    std::vector<T*> getVertices() {
         std::vector<T*> res;
         res.emplace_back(getFirstVertex());
         for(typename std::vector<E*>::iterator iter = edgesInOrder.begin(); iter != edgesInOrder.end(); iter++) {
             res.emplace_back(graph.getEdgeTarget(*iter));
         }
         return res;
     }

    T* getFirstVertex() {
         if(edgesInOrder.empty()) {
             return lastVertex;
         } else {
             return getGraph().getEdgeSource(edgesInOrder[0]);
         }
     }

     //返回指针需要delete
     uint8_t * getBases() {
         if(getEdges().empty()) {return graph.getAdditionalSequence(lastVertex);}

         uint8_t * res;
         uint8_t * bases = graph.getAdditionalSequence(graph.getEdgeSource(edgesInOrder[0]));
         res = bases;
         for(int i = 0; i < edgesInOrder.size(); i++) {
             uint8_t * tmp = new uint8_t[kmerSize+i+1];
             memcpy(tmp, res, kmerSize+i);
             delete[] bases;
             bases = graph.getAdditionalSequence(graph.getEdgeTarget(edgesInOrder[i]));
             memcpy(tmp+kmerSize+i, bases, 1);
             delete[] res;
             res = tmp;
         }
         delete[] bases;
         return res;
     }
};

#endif //MUTECT2CPP_MASTER_PATH_H
