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
    std::shared_ptr<V> source;
    std::shared_ptr<V> target;
    bool allowingMultipleEdges;


public:
    IntrusiveEdge(std::shared_ptr<V> source, std::shared_ptr<V> target) : source(source), target(target) {};
    std::shared_ptr<V> getSource() const {return source;}
    std::shared_ptr<V> getTarget() const {return target;}
};

template<class V, class E>
class DirectedSpecifics : public Specifics<V, E>{
private:
    static const long serialVersionUID = 8971725103718958232L;

    DirectedEdgeContainer<V, E> & getEdgeContainer(std::shared_ptr<V> vertex) {
        Mutect2Utils::validateArg(VertexSet.find(vertex) != VertexSet.end(), "no such vertex in graph");
        return vertexMapDirected.find(vertex)->second;
    }


    bool allowingLoops;
    bool allowingMultipleEdges;
    ArraySet<std::shared_ptr<V>> VertexSet;
    ArraySet<std::shared_ptr<E>> EdgeSet;

protected:
    std::map<std::shared_ptr<V>, DirectedEdgeContainer<V, E>> vertexMapDirected;

public:
    std::map<std::shared_ptr<E>, IntrusiveEdge<V>> edgeMap;

    DirectedSpecifics() = default;

    int degreeOf(std::shared_ptr<V> vertex) {throw std::invalid_argument("input argument");}

    void addVertex(std::shared_ptr<V> v) {
        if(v == nullptr) {throw std::invalid_argument("Null is not allowed here.");}
        if(containsVertex(v)) return;
        vertexMapDirected.insert(std::pair<std::shared_ptr<V>, DirectedEdgeContainer<V, E>>(v, DirectedEdgeContainer<V, E>()));
        VertexSet.insert(v);
    }

    ArraySet<std::shared_ptr<V>> & getVertexSet()  {
        return VertexSet;
    }

    ArraySet<std::shared_ptr<E>> getAllEdges(std::shared_ptr<V> sourceVertex, std::shared_ptr<V> targetVertex)  {
        ArraySet<std::shared_ptr<E>> edges;

        if(VertexSet.find(sourceVertex) != VertexSet.end() && VertexSet.find(targetVertex) != VertexSet.end()) {
            DirectedEdgeContainer<V, E> ec = getEdgeContainer(sourceVertex);
            typename std::vector<std::shared_ptr<E>>::iterator iter;
            for(iter = ec.outgoing.begin(); iter != ec.outgoing.end(); iter++){
                if(getEdgeTarget(*iter) == targetVertex)
                    edges.insert(*iter);
            }
        }
        return edges;
    }

    ArraySet<std::shared_ptr<V>> getAllTargets(std::shared_ptr<V> sourceVertex) {
        ArraySet<std::shared_ptr<V>> res;
        if(VertexSet.find(sourceVertex) != VertexSet.end()) {
            DirectedEdgeContainer<V, E> ec = getEdgeContainer(sourceVertex);
            typename std::vector<std::shared_ptr<E>>::iterator iter;
            for(iter = ec.outgoing.begin(); iter != ec.outgoing.end(); iter++){
                res.insert(getEdgeTarget(*iter));
            }
        }
        return res;
    }

    std::shared_ptr<V> getEdgeTarget(std::shared_ptr<E> e) {
        return edgeMap.find(e)->second.getTarget();
    }

    std::shared_ptr<E> getEdge(std::shared_ptr<V> sourceVertex, std::shared_ptr<V> targetVertex) {
        if(VertexSet.find(sourceVertex) != VertexSet.end() && VertexSet.find(targetVertex) != VertexSet.end()) {
            DirectedEdgeContainer<V, E> ec = getEdgeContainer(sourceVertex);
            typename std::vector<std::shared_ptr<E>>::iterator iter;
            for(iter = ec.outgoing.begin(); iter != ec.outgoing.end(); iter++){
                if(getEdgeTarget(*iter) == targetVertex)
                    return *iter;
            }
        }
        return nullptr;
    }

    void addEdgeToTouchingVertices(std::shared_ptr<E> e) {
        std::shared_ptr<V> source = getEdgeSource(e);
        std::shared_ptr<V> target = getEdgeTarget(e);

        getEdgeContainer(source).addOutgoingEdge(e);
        getEdgeContainer(target).addIncomingEdge(e);
    }

    std::shared_ptr<V> getEdgeSource(std::shared_ptr<E> e) {
        return edgeMap.find(e)->second.getSource();
    }

    ArraySet<std::shared_ptr<E>> edgesof(std::shared_ptr<V> vertex) {
        ArraySet<std::shared_ptr<E>> res = getEdgeContainer(vertex).incoming;
        ArraySet<std::shared_ptr<E>> tmp = getEdgeContainer(vertex).outgoing;

        typename std::vector<std::shared_ptr<E>>::iterator iter;
        for(iter = tmp.begin(); iter != tmp.end(); iter++) {
            res.insert(*iter);
        }
        if(allowingLoops) {
            ArraySet<std::shared_ptr<E>> loops = getAllEdges(vertex, vertex);
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

    int inDegreeOf(std::shared_ptr<V> vertex) {
        return getEdgeContainer(vertex).incoming.size();
    }

    int outDegreeOf(std::shared_ptr<V> vector) {
        return getEdgeContainer(vector).outgoing.size();
    }

    ArraySet<std::shared_ptr<E>> incomingEdgesOf(std::shared_ptr<V> vertex) {
        return getEdgeContainer(vertex).getUnmodifiableIncomingEdges();
    }


    ArraySet<std::shared_ptr<E>> outgoingEdgesOf(std::shared_ptr<V> vertex)  {
        return getEdgeContainer(vertex).getUnmodifiableOutgoingEdges();
    }

    void removeEdgeFromTouchingVertices(std::shared_ptr<E> e) {
        std::shared_ptr<V> source = getEdgeSource(e);
        std::shared_ptr<V> target = getEdgeTarget(e);

        getEdgeContainer(source).removeOutgoingEdge(e);
        getEdgeContainer(target).removeIncomingEdge(e);
    }

    bool addEdge(std::shared_ptr<V> sourceVertex, std::shared_ptr<V> targetVertex, std::shared_ptr<E> e) {
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
        edgeMap.insert(std::pair<std::shared_ptr<E>, IntrusiveEdge<V>>(e, intrusiveEdge));
        EdgeSet.insert(e);
        addEdgeToTouchingVertices(e);
        return true;
    }

    bool assertVertexExist(std::shared_ptr<V> v) {
        if(VertexSet.find(v) != VertexSet.end())
            return true;
        else
            throw std::invalid_argument("no such vertex in graph.");
    }

    bool containsVertex(std::shared_ptr<V> v) {
        return VertexSet.find(v) != VertexSet.end();
    }

    bool containsEdge(std::shared_ptr<E> e) {
        return EdgeSet.find(e) != EdgeSet.end();
    }

    bool isSource(std::shared_ptr<V> v) {
        Mutect2Utils::validateArg(v.get(), "Attempting to test a null vertex.");
        return inDegreeOf(v) == 0;
    }

    bool isSink(std::shared_ptr<V> v) {
        Mutect2Utils::validateArg(v.get(), "Attempting to test a null vertex.");
        return outDegreeOf(v) == 0;
    }

    std::shared_ptr<uint8_t[]> getAdditionalSequence(std::shared_ptr<V> v) {
        return v->getAdditionalSequence(isSource(v));
    }

    int getAdditionalSequenceLength(std::shared_ptr<V> v) {
        return v->getAdditionalSequenceLength(isSource(v));
    }

    bool removeEdge(std::shared_ptr<E> e) {
        if(containsEdge(e)) {
            removeEdgeFromTouchingVertices(e);
            edgeMap.erase(e);
            EdgeSet.erase(e);
            return true;
        } else {
            return false;
        }
    }

    std::shared_ptr<E> removeEdge(std::shared_ptr<V> sourceVertex, std::shared_ptr<V> targetVertex) {
        std::shared_ptr<E> e = getEdge(sourceVertex, targetVertex);

        if(e != nullptr) {
            removeEdgeFromTouchingVertices(e);
            EdgeSet.erase(e);
            edgeMap.erase(e);
        }
        return e;
    }

    std::shared_ptr<E> addEdge(std::shared_ptr<V> sourceVertex, std::shared_ptr<V> targetVertex) {
        assertVertexExist(sourceVertex);
        assertVertexExist(targetVertex);
        if(!allowingMultipleEdges
           && containsEdge(sourceVertex, targetVertex)) {
            return nullptr;
        }
        if (!allowingLoops && sourceVertex == targetVertex) {
            throw std::invalid_argument(LOOPS_NOT_ALLOWED);
        }
        std::shared_ptr<E> e = createEdge(sourceVertex, targetVertex);
        if(containsEdge(e)) {
            return nullptr;
        } else {
            edgeMap.insert(std::pair<std::shared_ptr<E>, IntrusiveEdge<V>>(e, IntrusiveEdge<V>(sourceVertex, targetVertex)));
            EdgeSet.insert(e);
            addEdgeToTouchingVertices(e);
            return e;
        }

    }

    bool containsEdge(std::shared_ptr<V> sourceVertex, std::shared_ptr<V> targetVertex) {
        return getEdge(sourceVertex, targetVertex) != nullptr;
    }

    virtual std::shared_ptr<E> createEdge(std::shared_ptr<V> sourceVertex, std::shared_ptr<V> targetVertex) {
        return std::shared_ptr<E>(new E());
    }

    virtual bool removeVertex(std::shared_ptr<V> v) {
        if(containsVertex(v)) {
            ArraySet<std::shared_ptr<E>> touchingEdgesList = edgesof(v);
            std::vector<std::shared_ptr<E>> edgesList;
            for(std::shared_ptr<E> e : touchingEdgesList)
                edgesList.template emplace_back(e);
            removeAllEdges(edgesList);
            vertexMapDirected.erase(v);
            VertexSet.erase(v);
            return true;
        } else {
            return false;
        }
    }

    bool removeAllEdges(std::vector<std::shared_ptr<E>> edges) {
        bool modified = false;
        for(auto e : edges) {
            modified |= removeEdge(e);
        }

        return modified;
    }

    bool removeAllEdges(std::list<std::shared_ptr<E>> edges) {
        bool modified = false;
        for(auto e : edges) {
            modified |= removeEdge(e);
        }

        return modified;
    }

    bool removeAllEdges(std::set<std::shared_ptr<E>> edges) {
        bool modified = false;
        for(auto e : edges) {
            modified |= removeEdge(e);
        }

        return modified;
    }

    bool removeAllVertices(std::vector<std::shared_ptr<V>> vertices) {
        bool modified = false;

        for(std::shared_ptr<V> v : vertices) {
            modified |= removeVertex(v);
        }
        return modified;
    }

    bool removeAllVertices(std::list<std::shared_ptr<V>> vertices) {
        bool modified = false;

        for(std::shared_ptr<V> v : vertices) {
            modified |= removeVertex(v);
        }
        return modified;
    }

    bool removeAllVertices(std::set<std::shared_ptr<V>> vertices) {
        bool modified = false;

        for(std::shared_ptr<V> v : vertices) {
            modified |= removeVertex(v);
        }
        return modified;
    }

    bool isRefSink(std::shared_ptr<V> v) {
        Mutect2Utils::validateArg(v.get() != nullptr, "Attempting to pull sequence from a null vertex.");

        for(std::shared_ptr<E> e : outgoingEdgesOf(v)){
            if(e->getIsRef())
                return false;
        }

        for(std::shared_ptr<E> e : incomingEdgesOf(v)){
            if(e->getIsRef())
                return true;
        }

        return VertexSet.size() == 1;
    }

    std::shared_ptr<E> incomingEdgeOf(std::shared_ptr<V> v) {
        Mutect2Utils::validateArg(v.get(), "Null is not allowed there");
        ArraySet<std::shared_ptr<E>> edgesSet = incomingEdgesOf(v);
        Mutect2Utils::validateArg(edgesSet.size() <= 1, "Cannot get a single incoming edge for a vertex with multiple incoming edges");
        return edgesSet.empty() ? nullptr : *edgesSet.begin();
    }

    std::shared_ptr<E> outgoingEdgeOf(std::shared_ptr<V> v) {
        Mutect2Utils::validateArg(v.get(), "Null is not allowed there");
        ArraySet<std::shared_ptr<E>> edgesSet = outgoingEdgesOf(v);
        Mutect2Utils::validateArg(edgesSet.size() <= 1, "Cannot get a single incoming edge for a vertex with multiple incoming edges");
        return edgesSet.empty() ? nullptr : *edgesSet.begin();
    }

    std::shared_ptr<V> getNextReferenceVertex(std::shared_ptr<V> v, bool allowNonRefPaths, std::shared_ptr<E> blacklistedEdge) {
        if(v == nullptr)
            return nullptr;
        ArraySet<std::shared_ptr<E>> outgoingEdges = outgoingEdgesOf(v);

        if(outgoingEdges.empty())
            return nullptr;

        for(std::shared_ptr<E> edgeToTest : outgoingEdges) {
            if(edgeToTest->getIsRef()) {
                return getEdgeTarget(edgeToTest);
            }
        }

        if(!allowNonRefPaths)
            return nullptr;

        std::vector<std::shared_ptr<E>> edges;
        for(std::shared_ptr<E> edgeToTest : outgoingEdges) {
            if(edgeToTest == blacklistedEdge) {
                edges.template emplace_back(edgeToTest);
            }
            if(edges.size() > 1)
                break;
        }
        return edges.size() == 1 ? getEdgeTarget(edges.at(0)) : nullptr;
    }

    std::shared_ptr<V> getPrevReferenceVertex(std::shared_ptr<V> v) {
        if(v == nullptr)
            return nullptr;
        ArraySet<std::shared_ptr<E>> edges = incomingEdgesOf(v);
        std::vector<std::shared_ptr<V>> allVertexs;
        for(std::shared_ptr<E> edge : edges) {
            std::shared_ptr<V> v = getEdgeSource(edge);
            if(isReferenceNode(v))
                allVertexs.template emplace_back(v);
        }
        return allVertexs.size() > 0 ? allVertexs.at(0) : nullptr;
    }

    bool isReferenceNode(std::shared_ptr<V> v) {
        Mutect2Utils::validateArg(v.get() != nullptr, "Attempting to test a null vertex.");
        for(std::shared_ptr<E> edge : edgesof(v)) {
            if(edge->getIsRef()) {
                return true;
            }
        }
        return VertexSet.size() == 1;
    }

    std::shared_ptr<V> getReferenceSourceVertex() {
        for(std::shared_ptr<V> vertex : VertexSet) {
            if(isRefSource(vertex))
                return vertex;
        }
        return nullptr;
    }

    std::shared_ptr<V> getReferenceSinkVertex() {
        for(std::shared_ptr<V> vertex : VertexSet) {
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
        ArraySet<std::shared_ptr<V>> onPathFromRefSource;
        BaseGraphIterator<V, E> sourceIter = BaseGraphIterator<V, E>(this, getReferenceSourceVertex(), false, true);
        while(sourceIter.hasNext()) {
            onPathFromRefSource.insert(sourceIter.next());
        }
        ArraySet<std::shared_ptr<V>> onPathFromRefSink;
        BaseGraphIterator<V, E> sinkIter = BaseGraphIterator<V, E>(this, getReferenceSinkVertex(), true, false);
        while(sinkIter.hasNext()) {
            onPathFromRefSink.insert(sinkIter.next());
        }
        //ArraySet<std::shared_ptr<V>> allVertex = getVertexSet();
        ArraySet<std::shared_ptr<V>> verticesToRemove;
        for(std::shared_ptr<V> v : VertexSet) {
            verticesToRemove.insert(v);
        }
        typename std::vector<std::shared_ptr<V>>::iterator iter = onPathFromRefSource.begin();
        for(; iter != onPathFromRefSource.end();) {
            if(onPathFromRefSink.find(*iter) == onPathFromRefSink.end()) {
                onPathFromRefSource.erase(iter);
            } else {
                iter++;
            }
        }
//        for(std::shared_ptr<V> v : onPathFromRefSource) {
//            if(onPathFromRefSink.find(v) == onPathFromRefSink.end()) {
//                onPathFromRefSource.erase(v);
//            }
//        }
        for(std::shared_ptr<V> v : onPathFromRefSource) {
            verticesToRemove.erase(v);
        }

        std::vector<std::shared_ptr<V>> vertices;
        for(std::shared_ptr<V> v : verticesToRemove) {
            vertices.emplace_back(v);
        }
        removeAllVertices(vertices);

        if ( getSinks().size() > 1 )
            throw std::length_error("Should have eliminated all but the reference sink");

        if ( getSources().size() > 1 )
            throw std::length_error("hould have eliminated all but the reference source");
    }

    ArraySet<std::shared_ptr<V>> incomingVerticesOf(std::shared_ptr<V> v) {
        Mutect2Utils::validateArg(v.get(), "Null is not allowed there.");
        ArraySet<std::shared_ptr<V>> ret;
        for(std::shared_ptr<E> e : incomingEdgesOf(v)) {
            ret.insert(getEdgeSource(e));
        }
        return ret;
    }

    ArraySet<std::shared_ptr<V>> outgoingVerticesOf(std::shared_ptr<V> v) {
        Mutect2Utils::validateArg(v.get(), "Null is not allowed there.");
        ArraySet<std::shared_ptr<V>> ret;
        for(std::shared_ptr<E> e : outgoingEdgesOf(v)) {
            ret.insert(getEdgeTarget(e));
        }
        return ret;
    }

    ArraySet<std::shared_ptr<V>> getSinks() {
        ArraySet<std::shared_ptr<V>> ret;
        for(std::shared_ptr<V> v : VertexSet) {
            if(isSink(v))
                ret.insert(v);
        }
        return ret;
    }

    ArraySet<std::shared_ptr<V>> getSources() {
        ArraySet<std::shared_ptr<V>> ret;
        for(std::shared_ptr<V> v : VertexSet) {
            if(isSource(v))
                ret.insert(v);
        }
        return ret;
    }

    void cleanNonRefPaths() {
        if(getReferenceSourceVertex() == nullptr || getReferenceSinkVertex() == nullptr ) {
            return;
        }
        std::set<std::shared_ptr<E>> edgesToCheck;
        for(std::shared_ptr<E> e : incomingEdgesOf(getReferenceSourceVertex())) {
            edgesToCheck.insert(e);
        }
        while(!edgesToCheck.empty()) {
            std::shared_ptr<E> e = *(edgesToCheck.begin());
            if(!e->getIsRef()) {
                for(std::shared_ptr<E> e : incomingEdgesOf(getEdgeSource(e))) {
                    edgesToCheck.insert(e);
                }
                removeEdge(e);
            }
            edgesToCheck.erase(e);
        }

        for(std::shared_ptr<E> e : outgoingEdgesOf(getReferenceSinkVertex())) {
            edgesToCheck.insert(e);
        }
        while(!edgesToCheck.empty()) {
            std::shared_ptr<E> e = *(edgesToCheck.begin());
            if(!e->getIsRef()) {
                for(std::shared_ptr<E> e : outgoingEdgesOf(getEdgeTarget(e))) {
                    edgesToCheck.insert(e);
                }
                removeEdge(e);
            }
            edgesToCheck.erase(e);
        }

        Specifics<V,E>::removeSingletonOrphanVertices();
    }


    bool isRefSource(std::shared_ptr<V> v){
        return Specifics<V,E>::isRefSource(v);
    }

    void removeVerticesNotConnectedToRefRegardlessOfEdgeDirection() {
        std::shared_ptr<V> refV = getReferenceSourceVertex();
        if(refV != nullptr) {
            BaseGraphIterator<V, E> iter = BaseGraphIterator<V, E>(this, refV, true, true);
            while(iter.hasNext()) {
                VertexSet.erase(iter.next());
            }
        }
        removeAllVertices(VertexSet.getArraySet());
    }

    void addOrUpdateEdge(std::shared_ptr<V> source, std::shared_ptr<V> target, std::shared_ptr<E> e) {
        Mutect2Utils::validateArg(source.get(), "source");
        Mutect2Utils::validateArg(target.get(), "target");
        Mutect2Utils::validateArg(e.get(), "edge");

        std::shared_ptr<E> prev = getEdge(source, target);
        if(prev != nullptr) {
            prev->add(*e);
        } else {
            addEdge(source, target, e);
        }
    }

    ArraySet<std::shared_ptr<E>> getEdgeSet() {
        return EdgeSet;
    }

    bool containsAllVertices(ArraySet<std::shared_ptr<V>> & vertices) {
        if(vertices.empty())
            throw std::invalid_argument("null vertex");
        for(std::shared_ptr<V> v : vertices) {
            if(!containsVertex(v))
                return false;
        }
        return true;
    }


};



#endif //MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H
