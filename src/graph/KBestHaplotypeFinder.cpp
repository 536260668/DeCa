//
// Created by 梦想家xixi on 2021/11/23.
//

#include "KBestHaplotypeFinder.h"
#include "BaseGraph/DFS_CycleDetect.h"
#include <queue>

KBestHaplotypeFinder::KBestHaplotypeFinder(SeqGraph *graph, ArraySet<SeqVertex *> & sources,
                                           ArraySet<SeqVertex *> & sinks) : graph(graph){
    Mutect2Utils::validateArg(graph, "graph cannot be null");
    Mutect2Utils::validateArg(!sources.empty(), "sources cannot be null");
    Mutect2Utils::validateArg(!sinks.empty(), "sinks cannot be null");
    Mutect2Utils::validateArg(graph->containsAllVertices(sources), "source does not belong to the graph");
    Mutect2Utils::validateArg(graph->containsAllVertices(sinks), "sink does not belong to the graph");

    this->graph =  DFS_CycleDetect<SeqVertex, BaseEdge>(*graph).detectCycles() ? removeCyclesAndVerticesThatDontLeadToSinks(graph, sources, sinks) : graph;
    this->sources = sources;
    this->sinks = sinks;
}

SeqGraph *
KBestHaplotypeFinder::removeCyclesAndVerticesThatDontLeadToSinks(SeqGraph *original, ArraySet<SeqVertex *> &sources,
                                                                 ArraySet<SeqVertex *> &sinks) {
    std::set<BaseEdge*> edgesToRemove;
    std::set<SeqVertex*> vertexToRemove;

    bool foundSomePath = false;
    for(SeqVertex* source : sources) {
        std::set<SeqVertex*> parentVertices;
        foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(original, source, sinks, edgesToRemove, vertexToRemove, parentVertices) || foundSomePath;
    }
    Mutect2Utils::validateArg(foundSomePath, "could not find any path from the source vertex to the sink vertex after removing cycles");
    Mutect2Utils::validateArg(!(edgesToRemove.empty() && vertexToRemove.empty()), "cannot find a way to remove the cycles");

    SeqGraph* result = original->clone();

    result->removeAllEdges(edgesToRemove);
    result->removeAllVertices(vertexToRemove);
    return result;
}

bool KBestHaplotypeFinder::findGuiltyVerticesAndEdgesToRemoveCycles(SeqGraph *graph, SeqVertex *currentVertex,
                                                                    ArraySet<SeqVertex *> &sinks,
                                                                    std::set<BaseEdge *> &edgesToRemove,
                                                                    std::set<SeqVertex *> &verticesToRemove,
                                                                    std::set<SeqVertex *> &parentVertices) {
    if(sinks.find(currentVertex) != sinks.end()) {
        return true;
    }
    ArraySet<BaseEdge*> outgoingEdges = graph->outgoingEdgesOf(currentVertex);
    parentVertices.insert(currentVertex);

    bool reachesSink = false;
    for(BaseEdge* edge : outgoingEdges) {
        SeqVertex* child = graph->getEdgeTarget(edge);
        if(parentVertices.find(child) != parentVertices.end()) {
            edgesToRemove.insert(edge);
        } else {
            bool childReachSink = findGuiltyVerticesAndEdgesToRemoveCycles(graph, child, sinks, edgesToRemove, verticesToRemove, parentVertices);
            reachesSink = reachesSink || childReachSink;
        }
    }
    if(!reachesSink) {
        verticesToRemove.insert(currentVertex);
    }
    return reachesSink;
}

KBestHaplotypeFinder::KBestHaplotypeFinder(SeqGraph *graph, SeqVertex *source, SeqVertex *sink) : graph(graph){
    sources.insert(source);
    sinks.insert(sink);
}

KBestHaplotypeFinder::KBestHaplotypeFinder(SeqGraph *graph) : graph(graph), sources(graph->getSources()), sinks(graph->getSinks()){}

std::vector<KBestHaplotype *> KBestHaplotypeFinder::findBestHaplotypes(int maxNumberOfHaplotypes) {
    std::vector<KBestHaplotype *> result;
    std::priority_queue<KBestHaplotype*, std::vector<KBestHaplotype*>, KBestHaplotypeComp> queue;
    for(SeqVertex* source : sources) {
        queue.push(new KBestHaplotype(source, *graph));
    }
    std::map<SeqVertex*, int> vertexCounts;
    for(SeqVertex* v : graph->getVertexSet()) {
        vertexCounts.insert(std::pair<SeqVertex*, int>(v, 0));
    }
    while(!queue.empty() && result.size() < maxNumberOfHaplotypes) {
        KBestHaplotype* pathToExtend = queue.top();
        queue.pop();
        SeqVertex* vertexToExtend = pathToExtend->getLastVertex();
        if(sinks.find(vertexToExtend) != sinks.end()) {
            result.emplace_back(pathToExtend);
        } else {
            ArraySet<BaseEdge*> outgoingEdges = graph->outgoingEdgesOf(vertexToExtend);
            int totalOutgoingMultiplicity = 0;
            for(BaseEdge* edge : outgoingEdges) {
                totalOutgoingMultiplicity += edge->getMultiplicity();
            }
            for(BaseEdge* edge : outgoingEdges) {
                SeqVertex* targetVertex = graph->getEdgeTarget(edge);
                int num = vertexCounts.at(targetVertex);
                num++;
                vertexCounts.find(targetVertex)->second = num;
                if(num < maxNumberOfHaplotypes) {
                    queue.push(new KBestHaplotype(pathToExtend, edge, totalOutgoingMultiplicity));
                }
            }
        }
    }
    return result;
}
