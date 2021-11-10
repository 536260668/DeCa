//
// Created by 梦想家xixi on 2021/10/21.
//

#ifndef MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H
#define MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H

#include "Specifics.h"
#include "Mutect2Utils.h"
#include "DirectedEdgeContainer.h"
#include <stdexcept>
#include <string>
#include <map>
#include <set>

static const std::string NOT_IN_DIRECTED_GRAPH = "no such operation in a directed graph";

static const std::string LOOPS_NOT_ALLOWED = "loops not allowed";

template<class V>
class IntrusiveEdge{
private:
    static const long serialVersionUID = 3258408452177932855L;
    V* source;
    V* target;
    bool allowingMultipleEdges;


public:
    IntrusiveEdge(V* source, V* target) : source(source), target(target) {};
    V* getSource() const {return source;}
    V* getTarget() const {return target;}
};

template<class V, class E>
class DirectedSpecifics : public Specifics<V, E>{
private:
    static const long serialVersionUID = 8971725103718958232L;

    DirectedEdgeContainer<V, E> & getEdgeContainer(V* vertex) {
        ArraySet<V*> allVector = this->getVertexSet();
        Mutect2Utils::validateArg(allVector.find(vertex) != allVector.end(), "no such vertex in graph");
        return vertexMapDirected.find(vertex)->second;
    }



    bool allowingLoops;
    bool allowingMultipleEdges;

protected:
    std::map<V*, DirectedEdgeContainer<V, E>> vertexMapDirected;

public:
    std::map<E*, IntrusiveEdge<V>> edgeMap;

    DirectedSpecifics() = default;

    int degreeOf(V* vertex) {throw std::invalid_argument("input argument");}

    void addVertex(V* v) {
        if(v == nullptr) {throw std::invalid_argument("Null is not allowed here.");}
        if(containsVertex(v)) return;
        vertexMapDirected.insert(std::pair<V*, DirectedEdgeContainer<V, E>>(v, DirectedEdgeContainer<V, E>()));
    }

    ArraySet<V*> getVertexSet()  {
        ArraySet<V*> res;
        typename  std::map<V*, DirectedEdgeContainer<V, E>>::iterator iter;
        for(iter = vertexMapDirected.begin(); iter != vertexMapDirected.end(); iter++) {
            res.insert(iter->first);
        }
        return res;
    }

    ArraySet<E*> getAllEdges(V* sourceVertex, V* targetVertex)  {
        ArraySet<E*> edges;
        ArraySet<V*> vertexs = getVertexSet();

        if(vertexs.find(sourceVertex) != vertexs.end() && vertexs.find(targetVertex) != vertexs.end()) {
            DirectedEdgeContainer<V, E> ec = getEdgeContainer(sourceVertex);
            typename std::vector<E*>::iterator iter;
            for(iter = ec.outgoing.begin(); iter != ec.outgoing.end(); iter++){
                if(getEdgeTarget(*iter) == targetVertex)
                    edges.insert(*iter);
            }
        }
        return edges;
    }

    ArraySet<V*> getAllTargets(V* sourceVertex) {
        ArraySet<V*> vertexs = getVertexSet();
        ArraySet<V*> res;
        if(vertexs.find(sourceVertex) != vertexs.end()) {
            DirectedEdgeContainer<V, E> ec = getEdgeContainer(sourceVertex);
            typename std::vector<E*>::iterator iter;
            for(iter = ec.outgoing.begin(); iter != ec.outgoing.end(); iter++){
                res.insert(getEdgeTarget(*iter));
            }
        }
        return res;
    }

    V* getEdgeTarget(E* e) {
        return edgeMap.find(e)->second.getTarget();
    }

    E* getEdge(V* sourceVertex, V* targetVertex) {
        ArraySet<V*> vertexs = getVertexSet();
        if(vertexs.find(sourceVertex) != vertexs.end() && vertexs.find(targetVertex) != vertexs.end()) {
            DirectedEdgeContainer<V, E> ec = getEdgeContainer(sourceVertex);
            typename std::vector<E*>::iterator iter;
            for(iter = ec.outgoing.begin(); iter != ec.outgoing.end(); iter++){
                if(getEdgeTarget(*iter) == targetVertex)
                    return *iter;
            }
        }
        return nullptr;
    }

    void addEdgeToTouchingVertices(E* e) {
        V* source = getEdgeSource(e);
        V* target = getEdgeTarget(e);

        getEdgeContainer(source).addOutgoingEdge(e);
        getEdgeContainer(target).addIncomingEdge(e);
    }

    V* getEdgeSource(E* e) {
        return edgeMap.find(e)->second.getSource();
    }

    ArraySet<E*> edgesof(V* vertex) {
        ArraySet<E*> res = getEdgeContainer(vertex).incoming;
        ArraySet<E*> tmp = getEdgeContainer(vertex).outgoing;

        typename std::vector<E*>::iterator iter;
        for(iter = tmp.begin(); iter != tmp.end(); iter++) {
            res.insert(*iter);
        }
        if(allowingLoops) {
            ArraySet<E*> loops = getAllEdges(vertex, vertex);
            for(iter = res.begin(); iter != res.end();) {
                if(loops.find(*iter) != loops.end()) {
                    loops.erase(iter);
                    res.erase(iter);
                }
                else
                    iter++;
            }
        }
        return res;
    }

    int inDegreeOf(V* vertex) {
        return getEdgeContainer(vertex).incoming.size();
    }

    int outDegreeOf(V* vector) {
        return getEdgeContainer(vector).outgoing.size();
    }

    ArraySet<E*> incomingEdgesOf(V* vertex) {
        return getEdgeContainer(vertex).getUnmodifiableIncomingEdges();
    }


    ArraySet<E*> outgoingEdgesOf(V* vertex)  {
        return getEdgeContainer(vertex).getUnmodifiableOutgoingEdges();
    }

    void removeEdgeFromTouchingVertices(E* e) {
        V* source = getEdgeSource(e);
        V* target = getEdgeTarget(e);

        getEdgeContainer(source).removeOutgoingEdge(e);
        getEdgeContainer(target).removeIncomingEdge(e);
    }

    bool addEdge(V* sourceVertex, V* targetVertex, E* e) {
        if(edgeMap.find(e) != edgeMap.end())
            return false;

        assertVertexExist(sourceVertex);
        assertVertexExist(targetVertex);

        if(!allowingMultipleEdges && edgeMap.find(getEdge(sourceVertex, targetVertex)) != edgeMap.end())
            return false;

        if(!allowingLoops && sourceVertex == targetVertex) {
            throw std::invalid_argument(LOOPS_NOT_ALLOWED);
        }

        IntrusiveEdge<V> intrusiveEdge(sourceVertex, targetVertex);
        edgeMap.insert(std::pair<E*, IntrusiveEdge<V>>(e, intrusiveEdge));
        addEdgeToTouchingVertices(e);
        return true;
    }

    bool assertVertexExist(V* v) {
        if(getVertexSet().find(v) != getVertexSet().end())
            return true;
        else
            throw std::invalid_argument("no such vertex in graph.");
    }

    bool containsVertex(V* v) {
        ArraySet<V*> res = getVertexSet();
        return res.find(v) != res.end();
    }

    bool containsEdge(E* e) {
        ArraySet<E*> res;
        for(typename std::map<E*, IntrusiveEdge<V>>::iterator iter = edgeMap.begin(); iter != edgeMap.end(); iter++)
            res.insert(iter->first);
        return res.find(e) != res.end();
    }

    bool isSource(V *v) {
        Mutect2Utils::validateArg(v, "Attempting to test a null vertex.");
        return inDegreeOf(v) == 0;
    }

    bool isSink(V *v) {
        Mutect2Utils::validateArg(v, "Attempting to test a null vertex.");
        return outDegreeOf(v) == 0;
    }

    uint8_t * getAdditionalSequence(V* v) {
        return v->getAdditionalSequence();
    }

    bool removeEdge(E* e) {
        if(containsEdge(e)) {
            removeEdgeFromTouchingVertices(e);
            edgeMap.erase(e);
            delete e;
            return true;
        } else {
            return false;
        }
    }

    virtual bool removeVertex(V *v) {
        if(containsVertex(v)) {
            ArraySet<E*> touchingEdgesList = edgesof(v);
            std::vector<E*> edgesList;
            for(E* e : touchingEdgesList)
                edgesList.template emplace_back(e);
            removeAllEdges(edgesList);
            vertexMapDirected.erase(v);
            delete v->getSequence();
            delete v;
            return true;
        } else {
            return false;
        }
    }

    bool removeAllEdges(std::vector<E*> edges) {
        bool modified = false;
        for(auto e : edges) {
            modified |= removeEdge(e);
        }

        return modified;
    }

    bool removeAllVertices(std::vector<V*> vertices) {
        bool modified = false;

        for(V *v : vertices) {
            modified |= removeVertex(v);
        }
        return modified;
    }

};



#endif //MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H
