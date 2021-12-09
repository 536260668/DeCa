//
// Created by 梦想家xixi on 2021/10/18.
//
#include "graph/Kmer.h"
#include "graph/MultiDeBruijnVertex.h"
#include "graph/ReadThreadingGraph.h"
#include "ReadThreadingGraph.h"
#include <iostream>
#include "path/AdaptiveChainPruner.h"
#include "graph/utils/GraphUtils.h"
#include "intel/smithwaterman/IntelSmithWaterman.h"
#include "graph/KBestHaplotypeFinder.h"
#include "read/CigarUtils.h"
#include "GenotypeLikelihoods.h"
#include "StringUtils.h"
//#include<string>


int main() {
    uint8_t * bases = new uint8_t[86]{65, 84, 65, 67, 65, 67, 67, 67, 71, 71, 67, 65, 67, 67, 67, 84, 71, 84, 67, 67, 84, 71, 71, 65, 67, 65, 67, 71, 67, 84, 71, 84, 84, 71, 71, 67, 67, 84, 71, 71, 65, 84, 67, 84, 71, 65, 71, 67, 67, 67, 84, 71, 71, 84, 71, 71, 65, 71, 71, 84, 67, 65, 65, 65, 71, 67, 67, 65, 67, 67, 84, 84, 84, 71, 71, 84, 84, 67, 84, 71, 67, 67, 65, 84, 84, 71};
    Kmer kmer(bases, 86);
    std::string key("abc");
    std::cout << (key == "abc") << std::endl;
//    Kmer subkmer = kmer.subKmer(10,10);
//    Kmer subkmer2 = kmer.subKmer(20,20);
//    int * diff = new int[20];
//    uint8_t * diff1 = new uint8_t[20];
//    int k = subkmer2.getDifferingPositions(subkmer, 10, diff, diff1);
//    Kmer subkmer3(subkmer2);
//    delete[] diff1;
//    delete[] diff;
//    delete[] bases;

    MultiDeBruijnVertex vertex(bases, 86, false);
    double i = -1.0/0.0;
    bool j = i > 1.0;
    smithwaterman_initial();
    CigarOperatorUtils::initial();
    std::cout << vertex.hasAmbiguousSequence() << std::endl;
    std::cout << vertex << std::endl;
    SequenceForKmers sequenceForKmers = {"chr1", bases, 0, 85, 1, false};
    ReadThreadingGraph graph = ReadThreadingGraph(10, 85, false, Kmer(nullptr, 0), 1);
    graph.setPending();
    graph.buildGraphIfNecessary();
    AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge> chainPruner = AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(0.001, 2.302585092994046, 100);
    chainPruner.pruneLowWeightChains(graph);
    graph.recoverDanglingTails(0, 4, false);
    graph.recoverDanglingHeads(0,4, false);

    graph.removePathsNotConnectedToRef();
    SeqGraph* seqGraph = graph.toSequenceGraph();
    seqGraph->cleanNonRefPaths();
    seqGraph->zipLinearChains();
    seqGraph->removeSingletonOrphanVertices();
    seqGraph->removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
    //GraphUtils::graphEquals(seqGraph, seqGraph);
    SeqVertex* source = seqGraph->getReferenceSourceVertex();
    SeqVertex* sink = seqGraph->getReferenceSinkVertex();
    KBestHaplotypeFinder finder = KBestHaplotypeFinder(seqGraph, source, sink);
    std::vector<KBestHaplotype*> res = finder.findBestHaplotypes(128);
    uint8_t * ref = new uint8_t[]{67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,65,67,67,67,84,65,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67};
    uint8_t * h = new uint8_t[]{67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,65,67,67,67,84,65,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67,67,84,65,65,67,67};
    Cigar* cigar = CigarUtils::calculateCigar(ref, 286, h, 285);
    delete[] bases;
    return 0;
}

