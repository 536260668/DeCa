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

class HaplotypeComp
{
public:
    bool operator()(const Haplotype* left, const Haplotype* right)
    {
        return (*left) < (*right);
    }
};

AssemblyResult *
ReadThreadingAssembler::getAssemblyResult(Haplotype *refHaplotype, int kmerSize, ReadThreadingGraph *rtgraph) {
    if(recoverDanglingBranches){
        rtgraph->recoverDanglingTails(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
        rtgraph->recoverDanglingHeads(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
    }

    if(removePathsNotConnectedToRef){
        rtgraph->removePathsNotConnectedToRef();
    }

    SeqGraph* initialSeqGraph = rtgraph->toSequenceGraph();
    if(justReturnRawGraph) {
        return new AssemblyResult(ASSEMBLED_SOME_VARIATION, initialSeqGraph, nullptr);
    }
    initialSeqGraph->cleanNonRefPaths();

    AssemblyResult* cleaned = cleanupSeqGraph(initialSeqGraph);
    Status status = cleaned->getStatus();
    AssemblyResult* ret = new AssemblyResult(status, cleaned->getGraph(), rtgraph);
    delete cleaned;
    return ret;
}

AssemblyResult *ReadThreadingAssembler::cleanupSeqGraph(SeqGraph *seqGraph) {
    seqGraph->zipLinearChains();
    seqGraph->removeSingletonOrphanVertices();
    seqGraph->removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
    seqGraph->simplifyGraph();
    if(seqGraph->getReferenceSinkVertex() == nullptr || seqGraph->getReferenceSourceVertex() == nullptr) {
        return new AssemblyResult(JUST_ASSEMBLED_REFERENCE, seqGraph, nullptr);
    }
    seqGraph->removePathsNotConnectedToRef();
    seqGraph->simplifyGraph();
    if(seqGraph->getVertexSet().size() == 1) {
        SeqVertex* complete = *(seqGraph->getVertexSet().begin());
        SeqVertex* dummy = new SeqVertex(nullptr, 0);
        seqGraph->addVertex(dummy);
        seqGraph->addEdge(complete, dummy, new BaseEdge(true, 0));
    }
    return new AssemblyResult(ASSEMBLED_SOME_VARIATION, seqGraph, nullptr);
}


std::vector<Haplotype *>
ReadThreadingAssembler::findBestPaths(const std::list<SeqGraph *>& graphs, Haplotype *refHaplotype, SimpleInterval *refLoc,
                                      SimpleInterval *activeRegionWindow,
                                      const std::map<SeqGraph *, AssemblyResult *>& assemblyResultByGraph, AssemblyResultSet* assemblyResultSet) const {
    std::set<Haplotype*, HaplotypeComp> returnHaplotypes;
    int activeRegionStart = refHaplotype->getAlignmentStartHapwrtRef();
    int failedCigars = 0;

    for(SeqGraph* graph : graphs) {
        SeqVertex* source = graph->getReferenceSourceVertex();
        SeqVertex* sink = graph->getReferenceSinkVertex();
        Mutect2Utils::validateArg(source != nullptr && sink != nullptr, "Both source and sink cannot be null");

        for(KBestHaplotype* kBestHaplotype : KBestHaplotypeFinder(graph, source, sink).findBestHaplotypes(numBestHaplotypesPerGraph)) {
            Haplotype* h = kBestHaplotype->getHaplotype();
            if(returnHaplotypes.find(h) == returnHaplotypes.end()) {
                if(kBestHaplotype->getIsReference()) {
                    refHaplotype->setScore(kBestHaplotype->getScore());
                }
                Cigar* cigar = CigarUtils::calculateCigar(refHaplotype->getBases(), refHaplotype->getLength(), h->getBases(), h->getLength());

                h->setCigar(cigar);
                h->setAlignmentStartHapwrtRef(activeRegionStart);
                h->setGenomeLocation(activeRegionWindow);
                returnHaplotypes.insert(h);
                assemblyResultSet->add(h, assemblyResultByGraph.at(graph));
            }
        }
    }
    if(returnHaplotypes.find(refHaplotype) == returnHaplotypes.end()) {
        returnHaplotypes.insert(refHaplotype);
    }
    //TODO:验证
    return {returnHaplotypes.begin(), returnHaplotypes.end()};
}

AssemblyResultSet *ReadThreadingAssembler::runLocalAssembly(AssemblyRegion *assemblyRegion, Haplotype *refHaplotype,
                                                            uint8_t *fullReferenceWithPadding, int refLength, SimpleInterval *refLoc,
                                                            ReadErrorCorrector *readErrorCorrector) {
    Mutect2Utils::validateArg(assemblyRegion, "Assembly engine cannot be used with a null AssemblyRegion.");
    Mutect2Utils::validateArg(refHaplotype, "Active region must have an extended location.");
    Mutect2Utils::validateArg(fullReferenceWithPadding, "fullReferenceWithPadding");
    Mutect2Utils::validateArg(refLoc, "refLoc");
    Mutect2Utils::validateArg(refLength == refLoc->size(), "Reference bases and reference loc must be the same size.");

    std::vector<SAMRecord> correctedReads;
    if(readErrorCorrector != nullptr) {
        //TODO::readErrorCorrector
        readErrorCorrector->addReadsToKmers(*assemblyRegion->getReads());
        correctedReads = *assemblyRegion->getReads();
    } else {
        correctedReads = *assemblyRegion->getReads();
    }
    std::vector<SeqGraph*> nonRefGraphs;
    AssemblyResultSet * resultSet = new AssemblyResultSet();
    resultSet->setRegionForGenotyping(assemblyRegion);
    resultSet->setFullReferenceWithPadding(fullReferenceWithPadding, refLength);
    resultSet->setPaddedReferenceLoc(refLoc);
    SimpleInterval activeRegionExtendedLocation = assemblyRegion->getExtendedSpan();
    refHaplotype->setGenomeLocation(&activeRegionExtendedLocation);
    resultSet->add(refHaplotype);
    std::map<SeqGraph*, AssemblyResult*> assemblyResultByGraph;
    for(AssemblyResult* result : assemble(correctedReads, refHaplotype)) {
        if(result->getStatus() == ASSEMBLED_SOME_VARIATION) {
            //TODO:do some QC on the graph
            assemblyResultByGraph.insert(std::pair<SeqGraph*, AssemblyResult*>(result->getGraph(), result));
            nonRefGraphs.emplace_back(result->getGraph());
        }
    }
    findBestPaths(nonRefGraphs, refHaplotype, refLoc, &activeRegionExtendedLocation, assemblyResultByGraph, resultSet);
    return resultSet;
}

std::vector<AssemblyResult *> ReadThreadingAssembler::assemble(std::vector<SAMRecord> &reads, Haplotype *refHaplotype) {
    std::vector<AssemblyResult *> results;
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

AssemblyResult *
ReadThreadingAssembler::createGraph(std::vector<SAMRecord> reads, Haplotype *refHaplotype, int kmerSize,
                                    bool allowLowComplexityGraphs, bool allowNonUniqueKmersInRef) {
    if(refHaplotype->getLength() < kmerSize) {
        return new AssemblyResult(FAILED, nullptr, nullptr);
    }
    SequenceForKmers tmp = {"ref", refHaplotype->getBases(), 0, refHaplotype->getLength(), 1, true};
    if(!allowNonUniqueKmersInRef && !ReadThreadingGraph::determineNonUniqueKmers(tmp, kmerSize).empty()) {
        return nullptr;
    }
    ReadThreadingGraph* rtgraph = new ReadThreadingGraph(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly, numPruningSamples);
    rtgraph->setThreadingStartOnlyAtExistingVertex(!recoverAllDanglingBranches);
    rtgraph->addSequence("ref", refHaplotype->getBases(), refHaplotype->getLength(), true);

    for(SAMRecord read : reads) {
        rtgraph->addRead(read);
    }
    //TODO:DELETE
    rtgraph->setPending();
    rtgraph->buildGraphIfNecessary();
    chainPruner->pruneLowWeightChains(*rtgraph);
    if(rtgraph->hasCycles()) {
        return nullptr;
    }

    if(! allowLowComplexityGraphs && rtgraph->isLowComplexity()) {
        return nullptr;
    }

    return getAssemblyResult(refHaplotype, kmerSize, rtgraph);
}

void ReadThreadingAssembler::addResult(std::vector<AssemblyResult *> &results, AssemblyResult *maybeNullResult) {
    if(maybeNullResult != nullptr) {
        results.emplace_back(maybeNullResult);
    }
}

int ReadThreadingAssembler::arrayMaxInt(std::vector<int> array) {
    Mutect2Utils::validateArg(!array.empty(), "Array size cannot be 0!");
    int ret = 0;
    for(int tmp : array) {
        if(ret < tmp)
            ret = tmp;
    }
    return ret;
}

std::vector<Haplotype *>
ReadThreadingAssembler::findBestPaths(const std::vector<SeqGraph *> &graphs, Haplotype *refHaplotype,
                                      SimpleInterval *refLoc, SimpleInterval *activeRegionWindow,
                                      const std::map<SeqGraph *, AssemblyResult *> &assemblyResultByGraph,
                                      AssemblyResultSet *assemblyResultSet) const {
    std::set<Haplotype*, HaplotypeComp> returnHaplotypes;
    int activeRegionStart = refHaplotype->getAlignmentStartHapwrtRef();
    int failedCigars = 0;

    for(SeqGraph* graph : graphs) {
        SeqVertex* source = graph->getReferenceSourceVertex();
        SeqVertex* sink = graph->getReferenceSinkVertex();
        Mutect2Utils::validateArg(source != nullptr && sink != nullptr, "Both source and sink cannot be null");

        for(KBestHaplotype* kBestHaplotype : KBestHaplotypeFinder(graph, source, sink).findBestHaplotypes(numBestHaplotypesPerGraph)) {
            Haplotype* h = kBestHaplotype->getHaplotype();
            if(returnHaplotypes.find(h) == returnHaplotypes.end()) {
                if(kBestHaplotype->getIsReference()) {
                    refHaplotype->setScore(kBestHaplotype->getScore());
                }
                Cigar* cigar = CigarUtils::calculateCigar(refHaplotype->getBases(), refHaplotype->getLength(), h->getBases(), h->getLength());

                h->setCigar(cigar);
                h->setAlignmentStartHapwrtRef(activeRegionStart);
                h->setGenomeLocation(activeRegionWindow);
                returnHaplotypes.insert(h);
                assemblyResultSet->add(h, assemblyResultByGraph.at(graph));
            }
        }
    }
    if(returnHaplotypes.find(refHaplotype) == returnHaplotypes.end()) {
        returnHaplotypes.insert(refHaplotype);
    }

    //TODO:验证
    return {returnHaplotypes.begin(), returnHaplotypes.end()};
}

ReadThreadingAssembler::ReadThreadingAssembler(int pruneFactor, int numPruningSamples, int numBestHaplotypesPerGraph,
                                               bool allowNonUniqueKmersInRef, std::vector<int> kmerSizes) : pruneFactor(pruneFactor), numPruningSamples(numPruningSamples), numBestHaplotypesPerGraph(numBestHaplotypesPerGraph)
                                               , allowNonUniqueKmersInRef(allowNonUniqueKmersInRef), kmerSizes(std::move(kmerSizes)){
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
