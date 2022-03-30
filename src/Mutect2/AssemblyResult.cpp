//
// Created by 梦想家xixi on 2021/11/15.
//

#include "AssemblyResult.h"
#include "Mutect2Utils.h"

AssemblyResult::AssemblyResult(Status status, const std::shared_ptr<SeqGraph>& graph, std::shared_ptr<ReadThreadingGraph> threadingGraph) : status(status), graph(graph), threadingGraph(threadingGraph){
    Mutect2Utils::validateArg(status == FAILED || graph != nullptr, "graph is null but status is not FAILED");
}

AssemblyResult::~AssemblyResult() {
	delete this->getThreadingGraph().get();
}
