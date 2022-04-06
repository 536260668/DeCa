//
// Created by 梦想家xixi on 2021/11/15.
//

#include "ReadThreadingAssembler.h"
#include <memory>
#include <utility>
#include "AssemblyResultSet.h"
#include "graph/KBestHaplotypeFinder.h"
#include "graph/ReadThreadingGraph.h"
#include "read/CigarUtils.h"
#include "AdaptiveChainPruner.h"
#include "boost/dynamic_bitset.hpp"


std::shared_ptr<AssemblyResult>
ReadThreadingAssembler::getAssemblyResult(std::shared_ptr<Haplotype> &refHaplotype, int kmerSize,
                                          const std::shared_ptr<ReadThreadingGraph> &rtgraph) const {
	//std::cout << rtgraph->getVertexSet().size() << rtgraph->getEdgeSet().size() << std::endl;
	if (recoverDanglingBranches) {
		rtgraph->recoverDanglingTails(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
		rtgraph->recoverDanglingHeads(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
	}
	//std::cout << rtgraph->getVertexSet().size() << rtgraph->getEdgeSet().size() << std::endl;
	if (removePathsNotConnectedToRef) {
		rtgraph->removePathsNotConnectedToRef();
	}

	std::shared_ptr<SeqGraph> initialSeqGraph = rtgraph->toSequenceGraph();
	if (justReturnRawGraph) {
		return std::make_shared<AssemblyResult>(ASSEMBLED_SOME_VARIATION, initialSeqGraph, nullptr);
	}
	initialSeqGraph->cleanNonRefPaths();

	std::shared_ptr<AssemblyResult> cleaned = cleanupSeqGraph(initialSeqGraph);
	Status status = cleaned->getStatus();
	std::shared_ptr<AssemblyResult> ret = std::make_shared<AssemblyResult>(status, cleaned->getGraph(), rtgraph);
	return ret;
}

std::shared_ptr<AssemblyResult> ReadThreadingAssembler::cleanupSeqGraph(const std::shared_ptr<SeqGraph> &seqGraph) {
	seqGraph->zipLinearChains();
	seqGraph->removeSingletonOrphanVertices();
	seqGraph->removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
	seqGraph->simplifyGraph();
	if (seqGraph->getReferenceSinkVertex() == nullptr || seqGraph->getReferenceSourceVertex() == nullptr) {
		return std::make_shared<AssemblyResult>(JUST_ASSEMBLED_REFERENCE, seqGraph, nullptr);
	}
	seqGraph->removePathsNotConnectedToRef();
	seqGraph->simplifyGraph();
	if (seqGraph->getVertexSet().size() == 1) {
		std::shared_ptr<SeqVertex> complete = *(seqGraph->getVertexSet().begin());
		std::shared_ptr<SeqVertex> dummy(new SeqVertex(nullptr, 0));
		seqGraph->addVertex(dummy);
		seqGraph->addEdge(complete, dummy, std::make_shared<BaseEdge>(true, 0));
	}
	return std::make_shared<AssemblyResult>(ASSEMBLED_SOME_VARIATION, seqGraph, nullptr);
}


std::vector<std::shared_ptr<Haplotype>>
ReadThreadingAssembler::findBestPaths(const std::list<std::shared_ptr<SeqGraph>> &graphs,
                                      std::shared_ptr<Haplotype> &refHaplotype,
                                      const std::shared_ptr<SimpleInterval> &refLoc,
                                      const std::shared_ptr<SimpleInterval> &activeRegionWindow,
                                      const std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>> &assemblyResultByGraph,
                                      std::shared_ptr<AssemblyResultSet> &assemblyResultSet) const {
	std::set<std::shared_ptr<Haplotype>, HaplotypeComp> returnHaplotypes;
	int activeRegionStart = refHaplotype->getAlignmentStartHapwrtRef();
	int failedCigars = 0;

	for (const std::shared_ptr<SeqGraph> &graph: graphs) {
		std::shared_ptr<SeqVertex> source = graph->getReferenceSourceVertex();
		std::shared_ptr<SeqVertex> sink = graph->getReferenceSinkVertex();
		Mutect2Utils::validateArg(source != nullptr && sink != nullptr, "Both source and sink cannot be null");

		for (const std::shared_ptr<KBestHaplotype> &kBestHaplotype: KBestHaplotypeFinder(graph, source,
		                                                                                 sink).findBestHaplotypes(
				numBestHaplotypesPerGraph)) {
			std::shared_ptr<Haplotype> h = kBestHaplotype->getHaplotype();
			if (returnHaplotypes.find(h) == returnHaplotypes.end()) {
				if (kBestHaplotype->getIsReference()) {
					refHaplotype->setScore(kBestHaplotype->getScore());
				}
				std::shared_ptr<Cigar> cigar = CigarUtils::calculateCigar(refHaplotype->getBases(),
				                                                          refHaplotype->getLength(), h->getBases(),
				                                                          h->getLength());

				h->setCigar(cigar);
				h->setAlignmentStartHapwrtRef(activeRegionStart);
				h->setGenomeLocation(activeRegionWindow);
				returnHaplotypes.insert(h);
				const std::shared_ptr<AssemblyResult> &tmp = assemblyResultByGraph.at(graph);
				assemblyResultSet->add(h, tmp);
			}
		}
	}
	if (returnHaplotypes.find(refHaplotype) == returnHaplotypes.end()) {
		returnHaplotypes.insert(refHaplotype);
	}
	//TODO:验证
	return {returnHaplotypes.begin(), returnHaplotypes.end()};
}

std::shared_ptr<AssemblyResultSet>
ReadThreadingAssembler::runLocalAssembly(const std::shared_ptr<AssemblyRegion> &assemblyRegion,
                                         std::shared_ptr<Haplotype> &refHaplotype,
                                         const std::shared_ptr<uint8_t[]> &fullReferenceWithPadding, int refLength,
                                         const std::shared_ptr<SimpleInterval> &refLoc,
                                         ReadErrorCorrector *readErrorCorrector) {
	Mutect2Utils::validateArg(assemblyRegion.get(), "Assembly engine cannot be used with a null AssemblyRegion.");
	Mutect2Utils::validateArg(refHaplotype.get(), "Active region must have an extended location.");
	Mutect2Utils::validateArg(fullReferenceWithPadding.get(), "fullReferenceWithPadding");
	Mutect2Utils::validateArg(refLoc.get(), "refLoc");
	Mutect2Utils::validateArg(refLength == refLoc->size(), "Reference bases and reference loc must be the same size.");

	std::vector<std::shared_ptr<SAMRecord>> correctedReads;
	if (readErrorCorrector != nullptr) {
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
	const std::shared_ptr<SimpleInterval> activeRegionExtendedLocation = assemblyRegion->getExtendedSpan();
	refHaplotype->setGenomeLocation(activeRegionExtendedLocation);
	resultSet->add(refHaplotype);
	std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>> assemblyResultByGraph;
	for (const std::shared_ptr<AssemblyResult> &result: assemble(correctedReads, refHaplotype)) {
		if (result->getStatus() == ASSEMBLED_SOME_VARIATION) {
			//TODO:do some QC on the graph
			assemblyResultByGraph.insert(std::make_pair(result->getGraph(), result));
			nonRefGraphs.emplace_back(result->getGraph());
		}
	}
	findBestPaths(nonRefGraphs, refHaplotype, refLoc, activeRegionExtendedLocation, assemblyResultByGraph, resultSet);
	return resultSet;
}

int
ReadThreadingAssembler::getMinKmerSize(std::shared_ptr<Haplotype> &refHaplotype, std::vector<int> candidateKmerSizes) {
	//std::string s = refHaplotype->getBaseString();
	std::shared_ptr<uint8_t[]> s = refHaplotype->getBases();
	int len = refHaplotype->getLength();
	int i, k = 0;

	//when candidateKmerSizes[k] <= 30, use Bit Operation (long long)
	uint8_t charToU8[len];
	for (i = 0; i < len; ++i) {
		uint8_t ch = s[i];
		if (ch == 'A') charToU8[i] = 0;
		else if (ch == 'C') charToU8[i] = 1;
		else if (ch == 'G') charToU8[i] = 2;
		else if (ch == 'T') charToU8[i] = 3;
		else throw std::invalid_argument("Found N in sequence when getting MinKmerSize!");
	}

	while (candidateKmerSizes[k] <= 30) {
		std::unordered_set<long long> valueSet;
		long long val = 0L, mask = (1L << (candidateKmerSizes[k] * 2)) - 1;
		for (i = 0; i < candidateKmerSizes[k]; ++i) val = (val << 2) | charToU8[i];
		valueSet.insert(val);
		for (i = candidateKmerSizes[k]; i < len; ++i) {
			val = ((val << 2) & mask) | charToU8[i];
			if (!valueSet.insert(val).second) {
				k++;
				break;
			}
		}
		if (i == len) return candidateKmerSizes[k];
	}

	//when candidateKmerSizes[k] > 30, use dynamic bitset
	uint8_t charToBitH[len], charToBitL[len];
	for (i = 0; i < len; ++i) {
		uint8_t ch = s[i];
		if (ch == 'A') charToBitH[i] = 0, charToBitL[i] = 0;
		else if (ch == 'C') charToBitH[i] = 0, charToBitL[i] = 1;
		else if (ch == 'G') charToBitH[i] = 1, charToBitL[i] = 0;
		else if (ch == 'T') charToBitH[i] = 1, charToBitL[i] = 1;
		else throw std::invalid_argument("Found N in sequence when getting MinKmerSize!");
	}

	int last_i = 0, j;
	while (k < candidateKmerSizes.size() - 1) {
		int bitSetSize = 2 * candidateKmerSizes[k];
		boost::dynamic_bitset<> s1(bitSetSize), s2(bitSetSize);
		for (j = last_i; j < last_i + candidateKmerSizes[k] - 1; ++j) {
			s1 <<= 2;
			s1[1] = charToBitH[j], s1[0] = charToBitL[j];
		}
		for (i = last_i + candidateKmerSizes[k] - 1; i < len; i++) {
			s1 <<= 2;
			s1[1] = charToBitH[i], s1[0] = charToBitL[i];
			//std::cout << "s1: " << s1 << std::endl;
			s2 = s1;
			for (j = i + 1; j < len; j++) { //s2[j] loop
				s2 <<= 2;
				s2[1] = charToBitH[j], s2[0] = charToBitL[j];
				//std::cout << "s2: " << s2 << std::endl;
				if (s1 == s2) { //match
					last_i = i - candidateKmerSizes[k] + 1;
					k++;
					if (last_i + candidateKmerSizes[k] >= len) return candidateKmerSizes[k];
					break;
				}
			}
			if (j != len) break; //match, solve next k
		}
		if (i == len) return candidateKmerSizes[k]; //not match, return
	}
	return candidateKmerSizes[k];
}

std::vector<std::shared_ptr<AssemblyResult>>
ReadThreadingAssembler::assemble(std::vector<std::shared_ptr<SAMRecord>> &reads,
                                 std::shared_ptr<Haplotype> &refHaplotype) {
	std::vector<std::shared_ptr<AssemblyResult>> results;
	if (kmerSizes.empty()) return results;

	std::sort(kmerSizes.begin(), kmerSizes.end());
	std::vector<int> tmpKmerSizes(kmerSizes);
	for (int numIterations = 0; numIterations < MAX_KMER_ITERATIONS_TO_ATTEMPT; ++numIterations) {
		tmpKmerSizes.push_back(*(tmpKmerSizes.end() - 1) + KMER_SIZE_ITERATION_INCREASE);
	}
	int minKmerSize = getMinKmerSize(refHaplotype, tmpKmerSizes);
	//std::cout << "minKmerSize = " << minKmerSize << std::endl;

	for (int kmerSize: kmerSizes) {
		if (kmerSize < minKmerSize && !allowNonUniqueKmersInRef) continue;
		addResult(results, createGraph(reads, refHaplotype, kmerSize, dontIncreaseKmerSizesForCycles));
	}

	if (results.empty() && !dontIncreaseKmerSizesForCycles) {
		int numIterations = 1, kmerSize = *(kmerSizes.end() - 1) + KMER_SIZE_ITERATION_INCREASE;
		while (numIterations < MAX_KMER_ITERATIONS_TO_ATTEMPT) {
			if (kmerSize >= minKmerSize || allowNonUniqueKmersInRef) {
				addResult(results, createGraph(reads, refHaplotype, kmerSize, false));
				if (!results.empty()) break;
			}
			numIterations++, kmerSize += KMER_SIZE_ITERATION_INCREASE;
		}
		if (numIterations == MAX_KMER_ITERATIONS_TO_ATTEMPT && results.empty())
			addResult(results, createGraph(reads, refHaplotype, kmerSize, true));
	}
	return results;
}

std::shared_ptr<AssemblyResult>
ReadThreadingAssembler::createGraph(const std::vector<std::shared_ptr<SAMRecord>> &reads,
                                    std::shared_ptr<Haplotype> &refHaplotype, int kmerSize,
                                    bool allowLowComplexityGraphs) {
	if (refHaplotype->getLength() < kmerSize) {
		return std::make_shared<AssemblyResult>(FAILED, nullptr, nullptr);
	}
	/* SequenceForKmers tmp = {"ref", refHaplotype->getBases(), 0, refHaplotype->getLength(), 1, true};
	 std::vector<std::shared_ptr<Kmer>>* res =ReadThreadingGraph::determineNonUniqueKmers(tmp, kmerSize);
	 if(!allowNonUniqueKmersInRef && !res->empty()) {
		 delete res;
		 return nullptr;
	 }
	 delete res;*/
	std::shared_ptr<ReadThreadingGraph> rtgraph = std::make_shared<ReadThreadingGraph>(kmerSize,
	                                                                                   debugGraphTransformations,
	                                                                                   minBaseQualityToUseInAssembly,
	                                                                                   numPruningSamples);
	rtgraph->setThreadingStartOnlyAtExistingVertex(!recoverDanglingBranches);
	rtgraph->addSequence("ref", refHaplotype->getBases(), refHaplotype->getLength(), true);
	rtgraph->reserveSpace(refHaplotype->getLength());

	for (std::shared_ptr<SAMRecord> read: reads) {
		rtgraph->addRead(read);
	}

	rtgraph->buildGraphIfNecessary();
	//std::cout << "1: " << rtgraph->getEdgeSet().size() << " " << rtgraph->getVertexSet().size() << std::endl;
	/*std::ofstream outfile1("./graph1.dot");
	outfile1 << "digraph G{" << std::endl;
	for (auto &v: rtgraph->getVertexSet()) {
		std::string s = reinterpret_cast<const char *>(v->getSequence().get());
		s = s.substr(0, v->getLength());
		//if (s==std::string("CCACAGCTCC")) std::cout<<"wdnmd ";
		outfile1 << "    " << s << ";" << std::endl;
	}
	for (auto &edge: rtgraph->edgeMap) {
		auto *a = edge.second.getSource().get();
		auto *b = edge.second.getTarget().get();
		std::string s1 = reinterpret_cast<const char *>(a->getSequence().get());
		std::string s2 = reinterpret_cast<const char *>(b->getSequence().get());
		s1 = s1.substr(0, a->getLength());
		s2 = s2.substr(0, b->getLength());
		outfile1 << "    " << s1 << " -> " << s2 << ";" << std::endl;
	}
	outfile1 << "}" << std::endl;
	outfile1.close();*/

	chainPruner->pruneLowWeightChains(rtgraph);
	//std::cout << "2: " << rtgraph->getEdgeSet().size() << " " << rtgraph->getVertexSet().size() << std::endl;
	/*std::ofstream outfile2("./graph2.dot");
		outfile2 << "digraph G{" << std::endl;
		for (auto &edge: rtgraph->edgeMap) {
			auto *a = edge.second.getSource().get();
			auto *b = edge.second.getTarget().get();
			std::string s1 = reinterpret_cast<const char *>(a->getSequence().get());
			std::string s2 = reinterpret_cast<const char *>(b->getSequence().get());
			s1 = s1.substr(0, a->getLength());
			s2 = s2.substr(0, b->getLength());
			outfile2 << "    " << s1 << " -> " << s2 << ";" << std::endl;
		}
		outfile2 << "}" << std::endl;
		outfile2.close();*/
	//    if(rtgraph->getVertexSet().size() == 292 && rtgraph->getEdgeSet().size() == 291)
	//        std::cout << " hello";
	//    std::ofstream outfile("/Users/bigdreamerxixi/data/1.txt", true);
	//    outfile << rtgraph->getVertexSet().size() << ", "<<rtgraph->getEdgeSet().size() << std::endl;
	//    outfile.close();

	if (rtgraph->hasCycles()) {
		//std::cout << kmerSize << " failed because hasCycles" << std::endl;
		return nullptr;
	}

	if (!allowLowComplexityGraphs && rtgraph->isLowComplexity()) {
		//std::cout << kmerSize << " failed because isLowComplexity" << std::endl;
		return nullptr;
	}
	//std::cout << kmerSize << std::endl;
	return getAssemblyResult(refHaplotype, kmerSize, rtgraph);
}

void ReadThreadingAssembler::addResult(std::vector<std::shared_ptr<AssemblyResult>> &results,
                                       const std::shared_ptr<AssemblyResult> &maybeNullResult) {
	if (maybeNullResult != nullptr) {
		results.emplace_back(maybeNullResult);
	}
}

std::vector<std::shared_ptr<Haplotype>>
ReadThreadingAssembler::findBestPaths(const std::vector<std::shared_ptr<SeqGraph>> &graphs,
                                      std::shared_ptr<Haplotype> &refHaplotype,
                                      const std::shared_ptr<SimpleInterval> &refLoc,
                                      const std::shared_ptr<SimpleInterval> &activeRegionWindow,
                                      const std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>> &assemblyResultByGraph,
                                      std::shared_ptr<AssemblyResultSet> &assemblyResultSet) const {
	std::set<std::shared_ptr<Haplotype>, HaplotypeComp> returnHaplotypes;
	int activeRegionStart = refHaplotype->getAlignmentStartHapwrtRef();
	int failedCigars = 0;

	for (const std::shared_ptr<SeqGraph> &graph: graphs) {
		std::shared_ptr<SeqVertex> source = graph->getReferenceSourceVertex();
		std::shared_ptr<SeqVertex> sink = graph->getReferenceSinkVertex();
		Mutect2Utils::validateArg(source != nullptr && sink != nullptr, "Both source and sink cannot be null");

		for (const std::shared_ptr<KBestHaplotype> &kBestHaplotype: KBestHaplotypeFinder(graph, source,
		                                                                                 sink).findBestHaplotypes(
				numBestHaplotypesPerGraph)) {
			std::shared_ptr<Haplotype> h = kBestHaplotype->getHaplotype();
			if (returnHaplotypes.find(h) == returnHaplotypes.end()) {
				if (kBestHaplotype->getIsReference()) {
					refHaplotype->setScore(kBestHaplotype->getScore());
				}
				std::shared_ptr<Cigar> cigar = CigarUtils::calculateCigar(refHaplotype->getBases(),
				                                                          refHaplotype->getLength(), h->getBases(),
				                                                          h->getLength());

				h->setCigar(cigar);
				h->setAlignmentStartHapwrtRef(activeRegionStart);
				h->setGenomeLocation(activeRegionWindow);
				returnHaplotypes.insert(h);
				const std::shared_ptr<AssemblyResult> &tmp = assemblyResultByGraph.at(graph);
				assemblyResultSet->add(h, tmp);
			}
		}
	}
	if (returnHaplotypes.find(refHaplotype) == returnHaplotypes.end()) {
		returnHaplotypes.insert(refHaplotype);
	}

	//TODO:验证
	return {returnHaplotypes.begin(), returnHaplotypes.end()};
}

ReadThreadingAssembler::ReadThreadingAssembler(int pruneFactor, int numPruningSamples, int numBestHaplotypesPerGraph,
                                               bool dontIncreaseKmerSizesForCycles,
                                               bool allowNonUniqueKmersInRef, std::vector<int> kmerSizes) : pruneFactor(
		pruneFactor), numPruningSamples(numPruningSamples), numBestHaplotypesPerGraph(numBestHaplotypesPerGraph),
                                                                                                            dontIncreaseKmerSizesForCycles(
		                                                                                                            dontIncreaseKmerSizesForCycles),
                                                                                                            allowNonUniqueKmersInRef(
		                                                                                                            allowNonUniqueKmersInRef),
                                                                                                            kmerSizes(
		                                                                                                            std::move(
				                                                                                                            kmerSizes)) {
	chainPruner = new AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(0.001, 2.302585092994046, 100);
	setMinDanglingBranchLength(4);
}

void ReadThreadingAssembler::setMinDanglingBranchLength(int minDanglingBranchLength) {
	this->minDanglingBranchLength = minDanglingBranchLength;
}

ReadThreadingAssembler::ReadThreadingAssembler(int maxAllowedPathsForReadThreadingAssembler, std::vector<int> kmerSizes,
                                               bool dontIncreaseKmerSizesForCycles, bool allowNonUniqueKmersInRef,
                                               int numPruningSamples, int pruneFactor, bool useAdaptivePruning,
                                               double initialErrorRateForPruning, double pruningLogOddsThreshold,
                                               int maxUnprunedVariants) : kmerSizes(std::move(kmerSizes)),
                                                                          dontIncreaseKmerSizesForCycles(
		                                                                          dontIncreaseKmerSizesForCycles),
                                                                          allowNonUniqueKmersInRef(
		                                                                          allowNonUniqueKmersInRef),
                                                                          numPruningSamples(numPruningSamples),
                                                                          pruneFactor(pruneFactor),
                                                                          numBestHaplotypesPerGraph(
		                                                                          maxAllowedPathsForReadThreadingAssembler) {
	Mutect2Utils::validateArg(maxAllowedPathsForReadThreadingAssembler >= 1,
	                          "numBestHaplotypesPerGraph should be >= 1");
	chainPruner = new AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(initialErrorRateForPruning,
	                                                                            pruningLogOddsThreshold,
	                                                                            maxUnprunedVariants);
}

ReadThreadingAssembler::~ReadThreadingAssembler() {
	delete chainPruner;
}

void ReadThreadingAssembler::setRecoverDanglingBranches(bool recoverDanglingBranches) {
	this->recoverDanglingBranches = recoverDanglingBranches;
}

void ReadThreadingAssembler::setRecoverAllDanglingBranches(bool recoverAllDanglingBranches) {
	this->recoverAllDanglingBranches = recoverAllDanglingBranches;
	recoverDanglingBranches = true;
}
