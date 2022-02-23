//
// Created by 梦想家xixi on 2021/11/15.
//

#include "ReadThreadingAssembler.h"
#include <utility>
#include "AssemblyResultSet.h"
#include "graph/KBestHaplotypeFinder.h"
#include "graph/ReadThreadingGraph.h"
#include "read/CigarUtils.h"
#include "AdaptiveChainPruner.h"



std::shared_ptr<AssemblyResult>
ReadThreadingAssembler::getAssemblyResult(std::shared_ptr<Haplotype>& refHaplotype, int kmerSize, const std::shared_ptr<ReadThreadingGraph>& rtgraph) {
    if(recoverDanglingBranches){
        rtgraph->recoverDanglingTails(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
        rtgraph->recoverDanglingHeads(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
    }

    if(removePathsNotConnectedToRef){
        rtgraph->removePathsNotConnectedToRef();
    }

    std::shared_ptr<SeqGraph> initialSeqGraph = rtgraph->toSequenceGraph();
    if(justReturnRawGraph) {
        return std::shared_ptr<AssemblyResult>(new AssemblyResult(ASSEMBLED_SOME_VARIATION, initialSeqGraph, nullptr));
    }
    initialSeqGraph->cleanNonRefPaths();

    std::shared_ptr<AssemblyResult> cleaned = cleanupSeqGraph(initialSeqGraph);
    Status status = cleaned->getStatus();
    std::shared_ptr<AssemblyResult> ret(new AssemblyResult(status, cleaned->getGraph(), rtgraph));
    return ret;
}

std::shared_ptr<AssemblyResult> ReadThreadingAssembler::cleanupSeqGraph(const std::shared_ptr<SeqGraph>&seqGraph) {
    seqGraph->zipLinearChains();
    seqGraph->removeSingletonOrphanVertices();
    seqGraph->removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
    seqGraph->simplifyGraph();
    if(seqGraph->getReferenceSinkVertex() == nullptr || seqGraph->getReferenceSourceVertex() == nullptr) {
        return std::shared_ptr<AssemblyResult>(new AssemblyResult(JUST_ASSEMBLED_REFERENCE, seqGraph, nullptr));
    }
    seqGraph->removePathsNotConnectedToRef();
    seqGraph->simplifyGraph();
    if(seqGraph->getVertexSet().size() == 1) {
        std::shared_ptr<SeqVertex> complete = *(seqGraph->getVertexSet().begin());
        std::shared_ptr<SeqVertex> dummy(new SeqVertex(nullptr, 0));
        seqGraph->addVertex(dummy);
        seqGraph->addEdge(complete, dummy, std::shared_ptr<BaseEdge>(new BaseEdge(true, 0)));
    }
    return std::shared_ptr<AssemblyResult>(new AssemblyResult(ASSEMBLED_SOME_VARIATION, seqGraph, nullptr));
}


std::vector<std::shared_ptr<Haplotype>>
ReadThreadingAssembler::findBestPaths(const std::list<std::shared_ptr<SeqGraph>>& graphs, std::shared_ptr<Haplotype>& refHaplotype, SimpleInterval *refLoc,
                                      SimpleInterval *activeRegionWindow,
                                      const std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>>& assemblyResultByGraph, std::shared_ptr<AssemblyResultSet>& assemblyResultSet) const {
    std::set<std::shared_ptr<Haplotype>, HaplotypeComp> returnHaplotypes;
    int activeRegionStart = refHaplotype->getAlignmentStartHapwrtRef();
    int failedCigars = 0;

    for(const std::shared_ptr<SeqGraph>& graph : graphs) {
        std::shared_ptr<SeqVertex> source = graph->getReferenceSourceVertex();
        std::shared_ptr<SeqVertex> sink = graph->getReferenceSinkVertex();
        Mutect2Utils::validateArg(source != nullptr && sink != nullptr, "Both source and sink cannot be null");

        for(const std::shared_ptr<KBestHaplotype>& kBestHaplotype : KBestHaplotypeFinder(graph, source, sink).findBestHaplotypes(numBestHaplotypesPerGraph)) {
            std::shared_ptr<Haplotype> h = kBestHaplotype->getHaplotype();
            if(returnHaplotypes.find(h) == returnHaplotypes.end()) {
                if(kBestHaplotype->getIsReference()) {
                    refHaplotype->setScore(kBestHaplotype->getScore());
                }
                std::shared_ptr<Cigar> cigar = CigarUtils::calculateCigar(refHaplotype->getBases(), refHaplotype->getLength(), h->getBases(), h->getLength());

                h->setCigar(cigar);
                h->setAlignmentStartHapwrtRef(activeRegionStart);
                h->setGenomeLocation(activeRegionWindow);
                returnHaplotypes.insert(h);
                std::shared_ptr<AssemblyResult> tmp = assemblyResultByGraph.at(graph);
                assemblyResultSet->add(h, tmp);
            }
        }
    }
    if(returnHaplotypes.find(refHaplotype) == returnHaplotypes.end()) {
        returnHaplotypes.insert(refHaplotype);
    }
    //TODO:验证
    return {returnHaplotypes.begin(), returnHaplotypes.end()};
}

std::shared_ptr<AssemblyResultSet> ReadThreadingAssembler::runLocalAssembly(const std::shared_ptr<AssemblyRegion>& assemblyRegion, std::shared_ptr<Haplotype> &refHaplotype,
                                                                            const std::shared_ptr<uint8_t[]>& fullReferenceWithPadding, int refLength, SimpleInterval *refLoc,
                                                            ReadErrorCorrector *readErrorCorrector) {
    Mutect2Utils::validateArg(assemblyRegion.get(), "Assembly engine cannot be used with a null AssemblyRegion.");
    Mutect2Utils::validateArg(refHaplotype.get(), "Active region must have an extended location.");
    Mutect2Utils::validateArg(fullReferenceWithPadding.get(), "fullReferenceWithPadding");
    Mutect2Utils::validateArg(refLoc, "refLoc");
    Mutect2Utils::validateArg(refLength == refLoc->size(), "Reference bases and reference loc must be the same size.");

    std::vector<std::shared_ptr<SAMRecord>> correctedReads;
    if(readErrorCorrector != nullptr) {
        //TODO::readErrorCorrector
        readErrorCorrector->addReadsToKmers(assemblyRegion->getReads());
        correctedReads = assemblyRegion->getReads();
    } else {
        correctedReads = assemblyRegion->getReads();
    }
    std::vector<std::shared_ptr<SeqGraph>> nonRefGraphs;
    std::shared_ptr<AssemblyResultSet> resultSet(new AssemblyResultSet());
    resultSet->setRegionForGenotyping(assemblyRegion);
    resultSet->setFullReferenceWithPadding(fullReferenceWithPadding, refLength);
    resultSet->setPaddedReferenceLoc(refLoc);
    SimpleInterval activeRegionExtendedLocation = assemblyRegion->getExtendedSpan();
    refHaplotype->setGenomeLocation(&activeRegionExtendedLocation);
    resultSet->add(refHaplotype);
    std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>> assemblyResultByGraph;
    for(const std::shared_ptr<AssemblyResult>& result : assemble(correctedReads, refHaplotype)) {
        if(result->getStatus() == ASSEMBLED_SOME_VARIATION) {
            //TODO:do some QC on the graph
            assemblyResultByGraph.insert(std::pair<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>>(result->getGraph(), result));
            nonRefGraphs.emplace_back(result->getGraph());
        }
    }
    findBestPaths(nonRefGraphs, refHaplotype, refLoc, &activeRegionExtendedLocation, assemblyResultByGraph, resultSet);
    return resultSet;
}

std::vector<std::shared_ptr<AssemblyResult>> ReadThreadingAssembler::assemble(std::vector<std::shared_ptr<SAMRecord>> &reads, std::shared_ptr<Haplotype> & refHaplotype) {
    std::vector<std::shared_ptr<AssemblyResult>> results;
    for(int kmerSize : kmerSizes) {
        addResult(results, createGraph(reads, refHaplotype, kmerSize, dontIncreaseKmerSizesForCycles, allowNonUniqueKmersInRef));
    }

    if(results.empty() && !dontIncreaseKmerSizesForCycles) {
        int kmerSize = arrayMaxInt(kmerSizes) + KMER_SIZE_ITERATION_INCREASE;
        int numIterations = 1;
        while(results.empty() && numIterations <= MAX_KMER_ITERATIONS_TO_ATTEMPT) {
            bool lastAttempt = numIterations == MAX_KMER_ITERATIONS_TO_ATTEMPT;
            addResult(results, createGraph(reads, refHaplotype, kmerSize, lastAttempt, lastAttempt));
            kmerSize += KMER_SIZE_ITERATION_INCREASE;
            numIterations++;
        }
    }
    return results;
}

std::shared_ptr<AssemblyResult>
ReadThreadingAssembler::createGraph(const std::vector<std::shared_ptr<SAMRecord>>& reads, std::shared_ptr<Haplotype>& refHaplotype, int kmerSize,
                                    bool allowLowComplexityGraphs, bool allowNonUniqueKmersInRef) {
    if(refHaplotype->getLength() < kmerSize) {
        return std::shared_ptr<AssemblyResult>(new AssemblyResult(FAILED, nullptr, nullptr));
    }
    SequenceForKmers tmp = {"ref", refHaplotype->getBases(), 0, refHaplotype->getLength(), 1, true};
    if(!allowNonUniqueKmersInRef && !ReadThreadingGraph::determineNonUniqueKmers(tmp, kmerSize).empty()) {
        return nullptr;
    }
    std::shared_ptr<ReadThreadingGraph> rtgraph(new ReadThreadingGraph(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, numPruningSamples));
    rtgraph->setThreadingStartOnlyAtExistingVertex(!recoverDanglingBranches);
    rtgraph->addSequence("ref", refHaplotype->getBases(), refHaplotype->getLength(), true);

    for(std::shared_ptr<SAMRecord> read : reads) {
        rtgraph->addRead(read);
    }

    rtgraph->buildGraphIfNecessary();


    chainPruner->pruneLowWeightChains(rtgraph);
//    if(rtgraph->getVertexSet().size() == 292 && rtgraph->getEdgeSet().size() == 291)
//        std::cout << " hello";
//    std::ofstream outfile("/Users/bigdreamerxixi/data/1.txt", true);
//    outfile << rtgraph->getVertexSet().size() << ", "<<rtgraph->getEdgeSet().size() << std::endl;
//    outfile.close();

    if(rtgraph->hasCycles()) {
        return nullptr;
    }

    if(! allowLowComplexityGraphs && rtgraph->isLowComplexity()) {
        return nullptr;
    }

    return getAssemblyResult(refHaplotype, kmerSize, rtgraph);
}

void ReadThreadingAssembler::addResult(std::vector<std::shared_ptr<AssemblyResult>> &results, const std::shared_ptr<AssemblyResult>& maybeNullResult) {
    if(maybeNullResult != nullptr) {
        results.emplace_back(maybeNullResult);
    }
}

int ReadThreadingAssembler::arrayMaxInt(const std::vector<int>& array) {
    Mutect2Utils::validateArg(!array.empty(), "Array size cannot be 0!");
    int ret = 0;
    for(int tmp : array) {
        if(ret < tmp)
            ret = tmp;
    }
    return ret;
}

std::vector<std::shared_ptr<Haplotype>>
ReadThreadingAssembler::findBestPaths(const std::vector<std::shared_ptr<SeqGraph>> &graphs, std::shared_ptr<Haplotype>& refHaplotype,
                                      SimpleInterval *refLoc, SimpleInterval *activeRegionWindow,
                                      const std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>> &assemblyResultByGraph,
                                      std::shared_ptr<AssemblyResultSet> &assemblyResultSet) const {
    std::set<std::shared_ptr<Haplotype>, HaplotypeComp> returnHaplotypes;
    int activeRegionStart = refHaplotype->getAlignmentStartHapwrtRef();
    int failedCigars = 0;

    for(const std::shared_ptr<SeqGraph>& graph : graphs) {
        std::shared_ptr<SeqVertex> source = graph->getReferenceSourceVertex();
        std::shared_ptr<SeqVertex> sink = graph->getReferenceSinkVertex();
        Mutect2Utils::validateArg(source != nullptr && sink != nullptr, "Both source and sink cannot be null");

        for(const std::shared_ptr<KBestHaplotype>& kBestHaplotype : KBestHaplotypeFinder(graph, source, sink).findBestHaplotypes(numBestHaplotypesPerGraph)) {
            std::shared_ptr<Haplotype> h = kBestHaplotype->getHaplotype();
            if(returnHaplotypes.find(h) == returnHaplotypes.end()) {
                if(kBestHaplotype->getIsReference()) {
                    refHaplotype->setScore(kBestHaplotype->getScore());
                }
                std::shared_ptr<Cigar> cigar = CigarUtils::calculateCigar(refHaplotype->getBases(), refHaplotype->getLength(), h->getBases(), h->getLength());

                h->setCigar(cigar);
                h->setAlignmentStartHapwrtRef(activeRegionStart);
                h->setGenomeLocation(activeRegionWindow);
                returnHaplotypes.insert(h);
                std::shared_ptr<AssemblyResult> tmp = assemblyResultByGraph.at(graph);
                assemblyResultSet->add(h, tmp);
            }
        }
    }
    if(returnHaplotypes.find(refHaplotype) == returnHaplotypes.end()) {
        returnHaplotypes.insert(refHaplotype);
    }

    //TODO:验证
    return {returnHaplotypes.begin(), returnHaplotypes.end()};
}

ReadThreadingAssembler::ReadThreadingAssembler(int pruneFactor, int numPruningSamples, int numBestHaplotypesPerGraph, bool dontIncreaseKmerSizesForCycles,
                                               bool allowNonUniqueKmersInRef, std::vector<int> kmerSizes) : pruneFactor(pruneFactor), numPruningSamples(numPruningSamples), numBestHaplotypesPerGraph(numBestHaplotypesPerGraph)
                                               , dontIncreaseKmerSizesForCycles(dontIncreaseKmerSizesForCycles),allowNonUniqueKmersInRef(allowNonUniqueKmersInRef), kmerSizes(std::move(kmerSizes)){
    chainPruner = new AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(0.001, 2.302585092994046, 100);
    setMinDanglingBranchLength(4);
}

void ReadThreadingAssembler::setMinDanglingBranchLength(int minDanglingBranchLength) {
    this->minDanglingBranchLength = minDanglingBranchLength;
}

ReadThreadingAssembler::ReadThreadingAssembler(int maxAllowedPathsForReadThreadingAssembler, std::vector<int> kmerSizes,
                                               bool dontIncreaseKmerSizesForCycles, bool allowNonUniqueKmersInRef,
                                               int numPruningSamples, int pruneFactor, bool useAdaptivePruning,
                                               double initialErrorRateForPruning, double pruningLogOddsThreshold, int maxUnprunedVariants) : kmerSizes(std::move(kmerSizes)), dontIncreaseKmerSizesForCycles(dontIncreaseKmerSizesForCycles), allowNonUniqueKmersInRef(allowNonUniqueKmersInRef),
                                               numPruningSamples(numPruningSamples), pruneFactor(pruneFactor), numBestHaplotypesPerGraph(maxAllowedPathsForReadThreadingAssembler){
    Mutect2Utils::validateArg(maxAllowedPathsForReadThreadingAssembler >= 1, "numBestHaplotypesPerGraph should be >= 1");
    chainPruner = new AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(initialErrorRateForPruning, pruningLogOddsThreshold, maxUnprunedVariants);
}

void ReadThreadingAssembler::setRecoverDanglingBranches(bool recoverDanglingBranches) {
    this->recoverDanglingBranches = recoverDanglingBranches;
}

void ReadThreadingAssembler::setRecoverAllDanglingBranches(bool recoverAllDanglingBranches) {
    this->recoverAllDanglingBranches = recoverAllDanglingBranches;
    recoverDanglingBranches = true;
}
