//
// Created by 梦想家xixi on 2021/11/15.
//

#include "SeqGraph.h"
#include <list>
#include <climits>
#include <cstring>
#include <memory>
#include "utils/MergeCommonSuffices.h"
#include "utils/MergeDiamonds.h"
#include "utils/MergeTails.h"
#include "utils/SplitCommonSuffices.h"
#include "utils/GraphUtils.h"

/*
std::shared_ptr<BaseEdge>
SeqGraph::createEdge(std::shared_ptr<SeqVertex> sourceVertex, std::shared_ptr<SeqVertex> targetVertrx) {
	return std::make_shared<BaseEdge>(false, 1);
}*/

bool SeqGraph::zipLinearChains() {
	std::vector<std::shared_ptr<SeqVertex>> zipStarts;
	for (auto &source: getVertexSet()) {
		if (isLinearChainStart(source)) {
			zipStarts.emplace_back(source);
		}
	}

	if (zipStarts.empty())
		return false;

	std::list<std::shared_ptr<SeqVertex>> linearChain;
	bool mergedOne = false;
	for (auto &zipStart: zipStarts) {
		linearChain = traceLinearChain(zipStart);
		mergedOne |= mergeLinearChain(linearChain);
	}
	return mergedOne;
}

bool SeqGraph::isLinearChainStart(const std::shared_ptr<SeqVertex> &source) {
	return outDegreeOf(source) == 1
	       && (inDegreeOf(source) != 1 || outDegreeOf(*(incomingVerticesOf(source).begin())) > 1);
}

std::list<std::shared_ptr<SeqVertex>> SeqGraph::traceLinearChain(const std::shared_ptr<SeqVertex> &zipStart) {
	std::list<std::shared_ptr<SeqVertex>> linearChain;
	linearChain.emplace_back(zipStart);

	bool lastIsRef = isReferenceNode(zipStart);
	std::shared_ptr<SeqVertex> last = zipStart;
	while (true) {
		if (outDegreeOf(last) != 1)
			break;

		std::shared_ptr<SeqVertex> target = getEdgeTarget(outgoingEdgeOf(last));
		if (inDegreeOf(target) != 1 || last == target)
			break;

		bool targetIsRef = isReferenceNode(target);
		if (lastIsRef != targetIsRef)
			break;
		linearChain.emplace_back(target);
		last = target;
		lastIsRef = targetIsRef;
	}
	return linearChain;
}

bool SeqGraph::mergeLinearChain(std::list<std::shared_ptr<SeqVertex>> &linearChain) {
	if (linearChain.empty())
		throw std::invalid_argument("BUG: cannot have linear chain with 0 elements");

	std::shared_ptr<SeqVertex> first = linearChain.front();
	std::shared_ptr<SeqVertex> last = linearChain.back();

	if (first == last)
		return false;

	std::shared_ptr<SeqVertex> addedVertex = mergeLinearChainVertices(linearChain);
	addVertex(addedVertex);

	for (auto &edge: outgoingEdgesOf(last)) {
		addEdge(addedVertex, getEdgeTarget(edge),
		        std::make_shared<BaseEdge>(edge->getIsRef(), edge->getMultiplicity()));
	}
	//std::unordered_set<std::shared_ptr<BaseEdge>> set1 = incomingEdgesOf(first);
	for (auto &edge: incomingEdgesOf(first)) {
		addEdge(getEdgeSource(edge), addedVertex,
		        std::make_shared<BaseEdge>(edge->getIsRef(), edge->getMultiplicity()));
	}
	removeAllVertices(linearChain);
	return true;
}


std::shared_ptr<SeqVertex> SeqGraph::mergeLinearChainVertices(std::list<std::shared_ptr<SeqVertex>> &vertices) {
	int length = 500;
	int start = 0;
	std::shared_ptr<uint8_t[]> tmp(new uint8_t[length]);
	for (auto &v: vertices) {
		int seqLength = v->getLength();
		while (start + seqLength >= length) {
			length <<= 1;
			std::shared_ptr<uint8_t[]> newtmp(new uint8_t[length]);
			memcpy(newtmp.get(), tmp.get(), start);
			tmp = newtmp;
		}
		memcpy(tmp.get() + start, v->getSequence().get(), seqLength);
		start += seqLength;
	}
	return std::make_shared<SeqVertex>(tmp, start);
}

void SeqGraph::simplifyGraph() {
	simplifyGraph(INT_MAX);
}

void SeqGraph::simplifyGraph(int maxCycles) {
	zipLinearChains();
	std::shared_ptr<SeqGraph> prevGraph = nullptr;
	for (int i = 0; i < maxCycles; i++) {
		if (i > MAX_REASONABLE_SIMPLIFICATION_CYCLES) {
			throw std::invalid_argument("Infinite loop detected in simplification routines for kmer graph");
		}
		if (!simplifyGraphOnce(i))
			break;
		if (i > 5) {
			if (prevGraph != nullptr && GraphUtils::graphEquals(prevGraph.get(), this))
				break;
		}
		prevGraph = std::make_shared<SeqGraph>(*this);
	}
}

bool SeqGraph::simplifyGraphOnce(int iteration) {
	bool didSomeWork = false;
	std::shared_ptr<SeqGraph> graph(new SeqGraph(*this));
	didSomeWork |= MergeDiamonds(graph).transformUntilComplete();
	didSomeWork |= MergeTails(graph).transformUntilComplete();
	didSomeWork |= SplitCommonSuffices(graph).transformUntilComplete();
	didSomeWork |= MergeCommonSuffices(graph).transformUntilComplete();
	didSomeWork |= zipLinearChains();
	return didSomeWork;
}

SeqGraph::SeqGraph(SeqGraph &seqGraph) : kmerSize(seqGraph.kmerSize), DirectedSpecifics<SeqVertex, BaseEdge>() {
	vertexMapDirected = seqGraph.vertexMapDirected;
	edgeMap = seqGraph.edgeMap;
}

std::shared_ptr<SeqGraph> SeqGraph::clone() {
	return std::make_shared<SeqGraph>(*this);
}

