//
// Created by 梦想家xixi on 2021/11/15.
//

#include "AssemblyResult.h"
#include "Mutect2Utils.h"

AssemblyResult::AssemblyResult(Status status, SeqGraph *graph, ReadThreadingGraph *threadingGraph) : status(status), graph(graph), threadingGraph(threadingGraph){
    Mutect2Utils::validateArg(status == FAILED || graph != nullptr, "graph is null but status is not FAILED");
}
