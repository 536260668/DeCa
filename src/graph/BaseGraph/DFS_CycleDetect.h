//
// Created by 梦想家xixi on 2021/11/5.
//

#ifndef MUTECT2CPP_MASTER_DFS_CYCLEDETECT_H
#define MUTECT2CPP_MASTER_DFS_CYCLEDETECT_H

#include "DirectedSpecifics.h"
#include <map>


template<class V, class E>
class DFS_CycleDetect {
public:
    bool detectCycles(){
        ArraySet<V*> allVertex = graph.getVertexSet();
        for(V* v : allVertex) {
            if(colorMap[v] == WHITE) {
                if(hasCycle(v))
                    return true;
            }
        }
        return false;
    }

    bool hasCycle(V* start){
        colorMap[start] = GREY;
        for(V* v : graph.getAllTargets(start)) {
            if(colorMap[v] == WHITE) {
                if(hasCycle(v))
                    return true;
            } else if(colorMap[v] == GREY) {
                return true;
            } else {continue;}
        }
        colorMap[start] = BLACK;
        return false;
    }

    DFS_CycleDetect(DirectedSpecifics<V, E> graph) : graph(graph) {
        ArraySet<V*> allVertex = graph.getVertexSet();
        for(V* v : allVertex) {
            colorMap.insert(std::pair<V*, color>(v, WHITE));
        }
    }

private:
    enum color{WHITE, GREY, BLACK};
    DirectedSpecifics<V, E> graph;
    std::map<V*, color> colorMap;
};


#endif //MUTECT2CPP_MASTER_DFS_CYCLEDETECT_H
