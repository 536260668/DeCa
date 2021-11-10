//
// Created by 梦想家xixi on 2021/10/21.
//

#ifndef MUTECT2CPP_MASTER_SPECIFICS_H
#define MUTECT2CPP_MASTER_SPECIFICS_H

#include "../set/ArraySet.h"
#include "Mutect2Utils.h"

template<class V, class E>
class Specifics {
private:
    static const long serialVersionUID = 785196247314761183L;

public:
    virtual void addVertex(V* vertex) = 0;

    virtual ArraySet<V*> getVertexSet() = 0;

    virtual ArraySet<E*> getAllEdges(V* sourceVertex, V* targetVertex) = 0;

    virtual E* getEdge(V* sourceVertex, V* targetVertex)  = 0;

    //Adds the specified edge to the edge containers of its source and target vertices.
    virtual void addEdgeToTouchingVertices(E* e) = 0;

    virtual int degreeOf(V* vertex)  = 0;

    virtual ArraySet<E*> edgesof(V* vertex) = 0;

    virtual int inDegreeOf(V* vertex)  = 0;

    virtual ArraySet<E*>  incomingEdgesOf(V* vertex) = 0;

    virtual int outDegreeOf(V* vertex)  = 0;

    virtual ArraySet<E*> outgoingEdgesOf(V* vertex) = 0;

    virtual void removeEdgeFromTouchingVertices(E* e) = 0;

    virtual bool isSource(V* v) = 0;

    virtual bool isSink(V* v) = 0;

    virtual bool isSingletonOrphan(V *v) {
        Mutect2Utils::validateArg(v, "v can not be null.");
        return inDegreeOf(v) == 0 && outDegreeOf(v) == 0 && !isRefSource(v);
    }

    virtual bool isRefSource(V *v) {
        Mutect2Utils::validateArg(v, "Attempting to pull sequence from a null vertex.");
        ArraySet<E*> incomingEdges = incomingEdgesOf(v);
        ArraySet<E*> outgoingEdges = outgoingEdgesOf(v);
        for(E *e : incomingEdges) {
            if(e->getIsRef())
                return false;
        }
        for(E *e : outgoingEdges) {
            if(e->getIsRef())
                return true;
        }

        return getVertexSet().size() == 1;
    }
    virtual void removeSingletonOrphanVertices() {
        std::vector<V*> toRemove;
        ArraySet<V*> allvertex = getVertexSet();
        typename std::vector<V*>::iterator viter;
        for(viter = allvertex.begin(); viter != allvertex.end(); viter++) {
            if(isSingletonOrphan(*viter)) {
                toRemove.template emplace_back(*viter);
            }
        }
        removeAllVertices(toRemove);
    }

    virtual bool removeAllVertices(std::vector<V*> vertices) = 0;

    virtual bool removeAllEdges(std::vector<E*> edges) = 0;
};


#endif //MUTECT2CPP_MASTER_SPECIFICS_H
