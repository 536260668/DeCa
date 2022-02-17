//
// Created by 梦想家xixi on 2021/11/23.
//

#include "KBestHaplotypeFinder.h"
#include "BaseGraph/DFS_CycleDetect.h"
#include <queue>

KBestHaplotypeFinder::KBestHaplotypeFinder(std::shared_ptr<SeqGraph> graph, ArraySet<std::shared_ptr<SeqVertex>> & sources,
                                           ArraySet<std::shared_ptr<SeqVertex>> & sinks) : graph(graph){
    Mutect2Utils::validateArg(graph.get(), "graph cannot be null");
    Mutect2Utils::validateArg(!sources.empty(), "sources cannot be null");
    Mutect2Utils::validateArg(!sinks.empty(), "sinks cannot be null");
    Mutect2Utils::validateArg(graph->containsAllVertices(sources), "source does not belong to the graph");
    Mutect2Utils::validateArg(graph->containsAllVertices(sinks), "sink does not belong to the graph");

    this->graph =  DFS_CycleDetect<SeqVertex, BaseEdge>(*graph).detectCycles() ? removeCyclesAndVerticesThatDontLeadToSinks(graph, sources, sinks) : graph;
    this->sources = sources;
    this->sinks = sinks;
}

std::shared_ptr<SeqGraph>
KBestHaplotypeFinder::removeCyclesAndVerticesThatDontLeadToSinks(std::shared_ptr<SeqGraph> original, ArraySet<std::shared_ptr<SeqVertex>> &sources,
                                                                 ArraySet<std::shared_ptr<SeqVertex>> &sinks) {
    std::set<std::shared_ptr<BaseEdge>> edgesToRemove;
    std::set<std::shared_ptr<SeqVertex>> vertexToRemove;

    bool foundSomePath = false;
    for(std::shared_ptr<SeqVertex> source : sources) {
        std::set<std::shared_ptr<SeqVertex>> parentVertices;
        foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(original, source, sinks, edgesToRemove, vertexToRemove, parentVertices) || foundSomePath;
    }
    Mutect2Utils::validateArg(foundSomePath, "could not find any path from the source vertex to the sink vertex after removing cycles");
    Mutect2Utils::validateArg(!(edgesToRemove.empty() && vertexToRemove.empty()), "cannot find a way to remove the cycles");

    std::shared_ptr<SeqGraph> result = original->clone();

    result->removeAllEdges(edgesToRemove);
    result->removeAllVertices(vertexToRemove);
    return result;
}

bool KBestHaplotypeFinder::findGuiltyVerticesAndEdgesToRemoveCycles(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex>currentVertex,
                                                                    ArraySet<std::shared_ptr<SeqVertex>> &sinks,
                                                                    std::set<std::shared_ptr<BaseEdge>> &edgesToRemove,
                                                                    std::set<std::shared_ptr<SeqVertex>> &verticesToRemove,
                                                                    std::set<std::shared_ptr<SeqVertex>> &parentVertices) {
    if(sinks.find(currentVertex) != sinks.end()) {
        return true;
    }
    ArraySet<std::shared_ptr<BaseEdge>> outgoingEdges = graph->outgoingEdgesOf(currentVertex);
    parentVertices.insert(currentVertex);

    bool reachesSink = false;
    for(std::shared_ptr<BaseEdge> edge : outgoingEdges) {
        std::shared_ptr<SeqVertex> child = graph->getEdgeTarget(edge);
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

KBestHaplotypeFinder::KBestHaplotypeFinder(std::shared_ptr<SeqGraph> graph, std::shared_ptr<SeqVertex>source, std::shared_ptr<SeqVertex>sink) : graph(graph){
    sources.insert(source);
    sinks.insert(sink);
}

KBestHaplotypeFinder::KBestHaplotypeFinder(std::shared_ptr<SeqGraph> graph) : graph(graph), sources(graph->getSources()), sinks(graph->getSinks()){}

std::vector<std::shared_ptr<KBestHaplotype>> KBestHaplotypeFinder::findBestHaplotypes(int maxNumberOfHaplotypes) {
    std::vector<std::shared_ptr<KBestHaplotype>> result;
    std::priority_queue<std::shared_ptr<KBestHaplotype>, std::vector<std::shared_ptr<KBestHaplotype>>, KBestHaplotypeComp> queue;
    for(std::shared_ptr<SeqVertex> source : sources) {
        queue.push(std::shared_ptr<KBestHaplotype>(new KBestHaplotype(source, graph)));
    }
    std::map<std::shared_ptr<SeqVertex>, int> vertexCounts;
    for(std::shared_ptr<SeqVertex> v : graph->getVertexSet()) {
        vertexCounts.insert(std::pair<std::shared_ptr<SeqVertex>, int>(v, 0));
    }
    while(!queue.empty() && result.size() < maxNumberOfHaplotypes) {
        std::shared_ptr<KBestHaplotype> pathToExtend = queue.top();
        queue.pop();
        std::shared_ptr<SeqVertex> vertexToExtend = pathToExtend->getLastVertex();
        if(sinks.find(vertexToExtend) != sinks.end()) {
            result.emplace_back(pathToExtend);
        } else {
            ArraySet<std::shared_ptr<BaseEdge>> outgoingEdges = graph->outgoingEdgesOf(vertexToExtend);
            int totalOutgoingMultiplicity = 0;
            for(std::shared_ptr<BaseEdge> edge : outgoingEdges) {
                totalOutgoingMultiplicity += edge->getMultiplicity();
            }
            for(std::shared_ptr<BaseEdge> edge : outgoingEdges) {
                std::shared_ptr<SeqVertex> targetVertex = graph->getEdgeTarget(edge);
                int num = vertexCounts.at(targetVertex);
                num++;
                vertexCounts.find(targetVertex)->second = num;
                if(num < maxNumberOfHaplotypes) {
                    queue.push(std::shared_ptr<KBestHaplotype>(new KBestHaplotype(pathToExtend, edge, totalOutgoingMultiplicity)));
                }
            }
        }
    }
    return result;
}
