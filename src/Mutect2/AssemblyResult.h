//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYRESULT_H
#define MUTECT2CPP_MASTER_ASSEMBLYRESULT_H

#include "ReadThreadingGraph.h"
#include "graph/SeqGraph.h"

enum Status{
    /** Something went wrong, and we couldn't produce a meaningful graph */
    FAILED,
    /** Assembly succeeded, but graph degenerated into just the reference sequence */
    JUST_ASSEMBLED_REFERENCE,
    /** Assembly succeeded, and the graph has some meaningful structure */
    ASSEMBLED_SOME_VARIATION
};

class AssemblyResult {
private:
    Status status;
    ReadThreadingGraph* threadingGraph;
    SeqGraph* graph;

public:
    AssemblyResult(Status status, SeqGraph* graph, ReadThreadingGraph* threadingGraph);

    ReadThreadingGraph* getThreadingGraph() {
        return threadingGraph;
    }

    Status getStatus() {
        return status;
    }

    SeqGraph* getGraph() {
        return graph;
    }

    int getKmerSize() {
        return graph->getKmerSize();
    }
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYRESULT_H
