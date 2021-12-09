//
// Created by 梦想家xixi on 2021/10/21.
//

#ifndef MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H
#define MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H

#include "Specifics.h"
#include "Mutect2Utils.h"
#include "DirectedEdgeContainer.h"
#include "BaseGraphIterator.h"
#include <stdexcept>
#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include "DFS_CycleDetect.h"

static const std::string NOT_IN_DIRECTED_GRAPH = "no such operation in a directed graph";

static const std::string LOOPS_NOT_ALLOWED = "loops not allowed";



template<class T, class E>
class BaseGraphIterator;


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
        return v->getAdditionalSequence(isSource(v));
    }

    int getAdditionalSequenceLength(V* v) {
        return v->getAdditionalSequenceLength(isSource(v));
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

    E* removeEdge(V* sourceVertex, V* targetVertex) {
        E* e = getEdge(sourceVertex, targetVertex);

        if(e != nullptr) {
            removeEdgeFromTouchingVertices(e);
            edgeMap.erase(e);
        }
        return e;
    }

    E* addEdge(V* sourceVertex, V* targetVertex) {
        assertVertexExist(sourceVertex);
        assertVertexExist(targetVertex);
        if(!allowingMultipleEdges
           && containsEdge(sourceVertex, targetVertex)) {
            return nullptr;
        }
        if (!allowingLoops && sourceVertex == targetVertex) {
            throw std::invalid_argument(LOOPS_NOT_ALLOWED);
        }
        E* e = createEdge(sourceVertex, targetVertex);
        if(containsEdge(e)) {
            return nullptr;
        } else {
            edgeMap.insert(std::pair<E*, IntrusiveEdge<V>>(e, IntrusiveEdge<V>(sourceVertex, targetVertex)));
            addEdgeToTouchingVertices(e);
            return e;
        }

    }

    bool containsEdge(V* sourceVertex, V* targetVertex) {
        return getEdge(sourceVertex, targetVertex) != nullptr;
    }

    virtual E* createEdge(V* sourceVertex, V* targetVertex) {
        return new E();
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

    bool removeAllEdges(std::list<E*> edges) {
        bool modified = false;
        for(auto e : edges) {
            modified |= removeEdge(e);
        }

        return modified;
    }

    bool removeAllEdges(std::set<E*> edges) {
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

    bool removeAllVertices(std::list<V*> vertices) {
        bool modified = false;

        for(V *v : vertices) {
            modified |= removeVertex(v);
        }
        return modified;
    }

    bool removeAllVertices(std::set<V*> vertices) {
        bool modified = false;

        for(V *v : vertices) {
            modified |= removeVertex(v);
        }
        return modified;
    }

    bool isRefSink(V* v) {
        Mutect2Utils::validateArg(v != nullptr, "Attempting to pull sequence from a null vertex.");

        for(E* e : outgoingEdgesOf(v)){
            if(e->getIsRef())
                return false;
        }

        for(E* e : incomingEdgesOf(v)){
            if(e->getIsRef())
                return true;
        }

        return getVertexSet().size() == 1;
    }

    E* incomingEdgeOf(V* v) {
        Mutect2Utils::validateArg(v, "Null is not allowed there");
        ArraySet<E*> edgesSet = incomingEdgesOf(v);
        Mutect2Utils::validateArg(edgesSet.size() <= 1, "Cannot get a single incoming edge for a vertex with multiple incoming edges");
        return edgesSet.empty() ? nullptr : *edgesSet.begin();
    }

    E* outgoingEdgeOf(V* v) {
        Mutect2Utils::validateArg(v, "Null is not allowed there");
        ArraySet<E*> edgesSet = outgoingEdgesOf(v);
        Mutect2Utils::validateArg(edgesSet.size() <= 1, "Cannot get a single incoming edge for a vertex with multiple incoming edges");
        return edgesSet.empty() ? nullptr : *edgesSet.begin();
    }

    V* getNextReferenceVertex(V* v, bool allowNonRefPaths, E* blacklistedEdge) {
        if(v == nullptr)
            return nullptr;
        ArraySet<E*> outgoingEdges = outgoingEdgesOf(v);

        if(outgoingEdges.empty())
            return nullptr;

        for(E* edgeToTest : outgoingEdges) {
            if(edgeToTest->getIsRef()) {
                return getEdgeTarget(edgeToTest);
            }
        }

        if(!allowNonRefPaths)
            return nullptr;

        std::vector<E*> edges;
        for(E* edgeToTest : outgoingEdges) {
            if(edgeToTest == blacklistedEdge) {
                edges.template emplace_back(edgeToTest);
            }
            if(edges.size() > 1)
                break;
        }
        return edges.size() == 1 ? getEdgeTarget(edges.at(0)) : nullptr;
    }

    V* getPrevReferenceVertex(V* v) {
        if(v == nullptr)
            return nullptr;
        ArraySet<E*> edges = incomingEdgesOf(v);
        std::vector<V*> allVertexs;
        for(E* edge : edges) {
            V* v = getEdgeSource(edge);
            if(isReferenceNode(v))
                allVertexs.template emplace_back(v);
        }
        return allVertexs.size() > 0 ? allVertexs.at(0) : nullptr;
    }

    bool isReferenceNode(V* v) {
        Mutect2Utils::validateArg(v != nullptr, "Attempting to test a null vertex.");
        for(E* edge : edgesof(v)) {
            if(edge->getIsRef()) {
                return true;
            }
        }
        return getVertexSet().size() == 1;
    }

    V* getReferenceSourceVertex() {
        for(V* vertex : getVertexSet()) {
            if(isRefSource(vertex))
                return vertex;
        }
        return nullptr;
    }

    V* getReferenceSinkVertex() {
        for(V* vertex : getVertexSet()) {
            if(isRefSink(vertex))
                return vertex;
        }
        return nullptr;
    }


    /**
     * Remove all vertices in the graph that aren't on a path from the reference source vertex to the reference sink vertex
     *
     * More aggressive reference pruning algorithm than removeVerticesNotConnectedToRefRegardlessOfEdgeDirection,
     * as it requires vertices to not only be connected by a series of directed edges but also prunes away
     * paths that do not also meet eventually with the reference sink vertex
     */
    void removePathsNotConnectedToRef() {
        if (getReferenceSourceVertex() == nullptr || getReferenceSinkVertex() == nullptr) {
            throw std::invalid_argument("Graph must have ref source and sink vertices");
        }
        ArraySet<V*> onPathFromRefSource;
        BaseGraphIterator<V, E> sourceIter = BaseGraphIterator<V, E>(this, getReferenceSourceVertex(), false, true);
        while(sourceIter.hasNext()) {
            onPathFromRefSource.insert(sourceIter.next());
        }
        ArraySet<V*> onPathFromRefSink;
        BaseGraphIterator<V, E> sinkIter = BaseGraphIterator<V, E>(this, getReferenceSinkVertex(), true, false);
        while(sinkIter.hasNext()) {
            onPathFromRefSink.insert(sinkIter.next());
        }
        ArraySet<V*> allVertex = getVertexSet();
        ArraySet<V*> verticesToRemove;
        for(V* v : allVertex) {
            verticesToRemove.insert(v);
        }
        for(V* v : onPathFromRefSource) {
            if(onPathFromRefSink.find(v) == onPathFromRefSink.end()) {
                onPathFromRefSource.erase(v);
            }
        }
        for(V* v : onPathFromRefSource) {
            verticesToRemove.erase(v);
        }

        std::vector<V*> vertices;
        for(V* v : verticesToRemove) {
            vertices.emplace_back(v);
        }
        removeAllVertices(vertices);

        if ( getSinks().size() > 1 )
            throw std::length_error("Should have eliminated all but the reference sink");

        if ( getSources().size() > 1 )
            throw std::length_error("hould have eliminated all but the reference source");
    }

    ArraySet<V*> incomingVerticesOf(V* v) {
        Mutect2Utils::validateArg(v, "Null is not allowed there.");
        ArraySet<V*> ret;
        for(E* e : incomingEdgesOf(v)) {
            ret.insert(getEdgeSource(e));
        }
        return ret;
    }

    ArraySet<V*> outgoingVerticesOf(V* v) {
        Mutect2Utils::validateArg(v, "Null is not allowed there.");
        ArraySet<V*> ret;
        for(E* e : outgoingEdgesOf(v)) {
            ret.insert(getEdgeTarget(e));
        }
        return ret;
    }

    ArraySet<V*> getSinks() {
        ArraySet<V*> ret;
        for(V* v : getVertexSet()) {
            if(isSink(v))
                ret.insert(v);
        }
        return ret;
    }

    ArraySet<V*> getSources() {
        ArraySet<V*> ret;
        for(V* v : getVertexSet()) {
            if(isSource(v))
                ret.insert(v);
        }
        return ret;
    }

    void cleanNonRefPaths() {
        if(getReferenceSourceVertex() == nullptr || getReferenceSinkVertex() == nullptr ) {
            return;
        }
        std::set<E*> edgesToCheck;
        for(E* e : incomingEdgesOf(getReferenceSourceVertex())) {
            edgesToCheck.insert(e);
        }
        while(!edgesToCheck.empty()) {
            E* e = *(edgesToCheck.begin());
            if(!e->getIsRef()) {
                for(E* e : incomingEdgesOf(getEdgeSource(e))) {
                    edgesToCheck.insert(e);
                }
                removeEdge(e);
            }
            edgesToCheck.erase(e);
        }

        for(E* e : outgoingEdgesOf(getReferenceSinkVertex())) {
            edgesToCheck.insert(e);
        }
        while(!edgesToCheck.empty()) {
            E* e = *(edgesToCheck.begin());
            if(!e->getIsRef()) {
                for(E* e : outgoingEdgesOf(getEdgeTarget(e))) {
                    edgesToCheck.insert(e);
                }
                removeEdge(e);
            }
            edgesToCheck.erase(e);
        }

        Specifics<V,E>::removeSingletonOrphanVertices();
    }


    bool isRefSource(V *v){
        return Specifics<V,E>::isRefSource(v);
    }

    void removeVerticesNotConnectedToRefRegardlessOfEdgeDirection() {
        ArraySet<V*> toRemove = getVertexSet();
        V* refV = getReferenceSourceVertex();
        if(refV != nullptr) {
            BaseGraphIterator<V, E> iter = BaseGraphIterator<V, E>(this, refV, true, true);
            while(iter.hasNext()) {
                toRemove.erase(iter.next());
            }
        }
        removeAllVertices(toRemove.getArraySet());
    }

    void addOrUpdateEdge(V* source, V* target, E* e) {
        Mutect2Utils::validateArg(source, "source");
        Mutect2Utils::validateArg(target, "target");
        Mutect2Utils::validateArg(e, "edge");

        E* prev = getEdge(source, target);
        if(prev != nullptr) {
            prev->add(*e);
        } else {
            addEdge(source, target, e);
        }
    }

    ArraySet<E*> getEdgeSet() {
        ArraySet<E*> ret;
        for(std::pair<E*, IntrusiveEdge<V>> edgePair : edgeMap) {
            ret.insert(edgePair.first);
        }
        return ret;
    }

    bool containsAllVertices(ArraySet<V*> & vertices) {
        if(vertices.empty())
            throw std::invalid_argument("null vertex");
        for(V* v : vertices) {
            if(!containsVertex(v))
                return false;
        }
        return true;
    }


};



#endif //MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H
