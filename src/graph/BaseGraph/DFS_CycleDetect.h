//
// Created by 梦想家xixi on 2021/11/5.
//

#ifndef MUTECT2CPP_MASTER_DFS_CYCLEDETECT_H
#define MUTECT2CPP_MASTER_DFS_CYCLEDETECT_H

#include "DirectedSpecifics.h"
#include <map>

template<class V, class E>
class DirectedSpecifics;

template<class V, class E>
class DFS_CycleDetect {
public:
    bool detectCycles(){
        std::unordered_set<std::shared_ptr<V>>& allVertex = graph->getVertexSet();
        for(const std::shared_ptr<V> & v : allVertex) {
            if(colorMap[v] == WHITE) {
                if(hasCycle(v))
                    return true;
            }
        }
        return false;
    }

    bool hasCycle(const std::shared_ptr<V> & start){
        colorMap[start] = GREY;
        for(const std::shared_ptr<V> & v : graph->getAllTargets(start)) {
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

    DFS_CycleDetect(DirectedSpecifics<V, E> * graph) : graph(graph) {
        std::unordered_set<std::shared_ptr<V>>& allVertex = graph->getVertexSet();
        for(const std::shared_ptr<V> & v : allVertex) {
            colorMap.insert(std::pair<std::shared_ptr<V>, color>(v, WHITE));
        }
    }

private:
    enum color{WHITE, GREY, BLACK};
    DirectedSpecifics<V, E> * graph;
    std::map<std::shared_ptr<V>, color> colorMap;
};


#endif //MUTECT2CPP_MASTER_DFS_CYCLEDETECT_H
