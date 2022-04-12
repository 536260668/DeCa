//
// Created by 梦想家xixi on 2021/10/21.
//

#ifndef MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H
#define MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H

#include "Specifics.h"
#include "Mutect2Utils.h"
#include "DirectedEdgeContainer.h"
#include "BaseGraphIterator.h"
#include <stdexcept>
#include <string>
#include <deque>
#include <list>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include "DFS_CycleDetect.h"

static const std::string NOT_IN_DIRECTED_GRAPH = "no such operation in a directed graph";

static const std::string LOOPS_NOT_ALLOWED = "loops not allowed";


template<class T, class E>
class BaseGraphIterator;


template<class V>
class IntrusiveEdge {
private:
	static const long serialVersionUID = 3258408452177932855L;
	std::shared_ptr<V> source;
	std::shared_ptr<V> target;
	bool allowingMultipleEdges;


public:
	IntrusiveEdge(std::shared_ptr<V> source, std::shared_ptr<V> target) : source(source), target(target) {};

	std::shared_ptr<V> getSource() const { return source; }

	std::shared_ptr<V> getTarget() const { return target; }
};

template<class V, class E>
class DirectedSpecifics : public Specifics<V, E> {
private:
	static const long serialVersionUID = 8971725103718958232L;

	DirectedEdgeContainer<E> &getEdgeContainer(const std::shared_ptr<V> &vertex) {
		auto miter = vertexMapDirected.find(vertex);
		if (miter == vertexMapDirected.end()) {
			throw std::invalid_argument("no such vertex in graph");
		}
		return miter->second;
	}

	bool allowingLoops{};
	bool allowingMultipleEdges{};
	std::unordered_set<std::shared_ptr<V>> VertexSet;
	std::unordered_set<std::shared_ptr<E>> EdgeSet;

protected:
	std::unordered_map<std::shared_ptr<V>, DirectedEdgeContainer<E>> vertexMapDirected;

public:
	std::unordered_map<std::shared_ptr<E>, IntrusiveEdge<V>> edgeMap;

	DirectedSpecifics() = default;

	~DirectedSpecifics() = default;

	int degreeOf(const std::shared_ptr<V> &vertex) { throw std::invalid_argument("input argument"); }

	void addVertex(const std::shared_ptr<V> &v) {
		if (v == nullptr)
			throw std::invalid_argument("Null is not allowed here.");
		if (containsVertex(v))
			return;
		vertexMapDirected.insert(std::make_pair(v, DirectedEdgeContainer<E>()));
		VertexSet.insert(v);
	}

	std::unordered_set<std::shared_ptr<V>> &getVertexSet() {
		return VertexSet;
	}

	std::unordered_set<std::shared_ptr<E>>
	getAllEdges(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		std::unordered_set<std::shared_ptr<E>> edges;

		if (VertexSet.find(sourceVertex) != VertexSet.end() && VertexSet.find(targetVertex) != VertexSet.end()) {
			const DirectedEdgeContainer<E> &ec = getEdgeContainer(sourceVertex);
			for (auto iter = ec.outgoing.begin(); iter != ec.outgoing.end(); iter++) {
				if (getEdgeTarget(*iter) == targetVertex)
					edges.insert(*iter);
			}
		}
		return edges;
	}

	std::unordered_set<std::shared_ptr<V>> getAllTargets(const std::shared_ptr<V> &sourceVertex) {
		std::unordered_set<std::shared_ptr<V>> res;
		if (VertexSet.find(sourceVertex) != VertexSet.end()) {
			const DirectedEdgeContainer<E> &ec = getEdgeContainer(sourceVertex);

			for (auto iter = ec.outgoing.begin(); iter != ec.outgoing.end(); iter++) {
				res.insert(getEdgeTarget(*iter));
			}
		}
		return res;
	}

	std::shared_ptr<V> getEdgeTarget(const std::shared_ptr<E> &e) {
		return edgeMap.find(e)->second.getTarget();
	}

	std::shared_ptr<E> getEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		if (VertexSet.find(sourceVertex) != VertexSet.end() && VertexSet.find(targetVertex) != VertexSet.end()) {
			const DirectedEdgeContainer<E> &ec = getEdgeContainer(sourceVertex);
			//typename std::unordered_set<std::shared_ptr<E>>::iterator iter;
			for (auto iter = ec.outgoing.begin(); iter != ec.outgoing.end(); iter++) {
				if (getEdgeTarget(*iter) == targetVertex)
					return *iter;
			}
		}
		return nullptr;
	}

	/*void addEdgeToTouchingVertices(const std::shared_ptr<E> &e) {
		std::shared_ptr<V> source = getEdgeSource(e);
		std::shared_ptr<V> target = getEdgeTarget(e);

		getEdgeContainer(source).addOutgoingEdge(e);
		getEdgeContainer(target).addIncomingEdge(e);
	}*/

	std::shared_ptr<V> getEdgeSource(std::shared_ptr<E> e) {
		return edgeMap.find(e)->second.getSource();
	}

	std::unordered_set<std::shared_ptr<E>> edgesof(const std::shared_ptr<V> &vertex) {
		std::unordered_set<std::shared_ptr<E>> res = getEdgeContainer(vertex).incoming;
		const std::unordered_set<std::shared_ptr<E>> &outgoing = getEdgeContainer(vertex).outgoing;
		res.reserve(res.size() + outgoing.size());
		for (const std::shared_ptr<E> &e: outgoing) {
			res.insert(e);
		}

		if (allowingLoops) {
			std::unordered_set<std::shared_ptr<E>> loops = getAllEdges(vertex, vertex);
			for (typename std::unordered_set<std::shared_ptr<E>>::iterator iter = res.begin();
			     iter != res.end(); iter++) {
				if (loops.find(*iter) != loops.end()) {
					loops.erase(iter);
					res.erase(*iter);
					//toRemove.insert(*iter);
				}
			}
		}
		return res;
	}

	int inDegreeOf(const std::shared_ptr<V> &vertex) {
		return getEdgeContainer(vertex).incoming.size();
	}

	int outDegreeOf(const std::shared_ptr<V> &vector) {
		return getEdgeContainer(vector).outgoing.size();
	}

	std::unordered_set<std::shared_ptr<E>> &incomingEdgesOf(const std::shared_ptr<V> &vertex) {
		return getEdgeContainer(vertex).getUnmodifiableIncomingEdges();
	}


	std::unordered_set<std::shared_ptr<E>> &outgoingEdgesOf(const std::shared_ptr<V> &vertex) {
		return getEdgeContainer(vertex).getUnmodifiableOutgoingEdges();
	}

	void removeEdgeFromTouchingVertices(const std::shared_ptr<E> &e) {
		std::shared_ptr<V> source = getEdgeSource(e);
		std::shared_ptr<V> target = getEdgeTarget(e);

		getEdgeContainer(source).removeOutgoingEdge(e);
		getEdgeContainer(target).removeIncomingEdge(e);
	}

	bool addEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex,
	             const std::shared_ptr<E> &e) {
		if (containsEdge(e) || !containsVertex(sourceVertex) || !containsVertex(targetVertex))
			return false;

		if (!allowingMultipleEdges && edgeMap.find(getEdge(sourceVertex, targetVertex)) != edgeMap.end())
			return false;

		if (!allowingLoops && sourceVertex == targetVertex) {
			throw std::invalid_argument(LOOPS_NOT_ALLOWED);
		}

		edgeMap.insert(std::make_pair(e, IntrusiveEdge<V>(sourceVertex, targetVertex)));
		EdgeSet.insert(e);
		getEdgeContainer(sourceVertex).addOutgoingEdge(e);
		getEdgeContainer(targetVertex).addIncomingEdge(e);
		return true;
	}

	bool assertVertexExist(const std::shared_ptr<V> &v) {
		if (vertexMapDirected.find(v) != vertexMapDirected.end())
			return true;
		throw std::invalid_argument("no such vertex in graph.");
	}

	bool containsVertex(const std::shared_ptr<V> &v) {
		return vertexMapDirected.find(v) != vertexMapDirected.end();
	}

	bool containsEdge(const std::shared_ptr<E> &e) {
		return edgeMap.find(e) != edgeMap.end();
	}

	bool isSource(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		return inDegreeOf(v) == 0;
	}

	bool isSink(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		return outDegreeOf(v) == 0;
	}

	std::shared_ptr<uint8_t[]> getAdditionalSequence(const std::shared_ptr<V> &v) {
		return v->getAdditionalSequence(isSource(v));
	}

	int getAdditionalSequenceLength(const std::shared_ptr<V> &v) {
		return v->getAdditionalSequenceLength(isSource(v));
	}

	bool removeEdge(const std::shared_ptr<E> &e) {
		if (containsEdge(e)) {
			removeEdgeFromTouchingVertices(e);
			edgeMap.erase(e);
			EdgeSet.erase(e);
			return true;
		}
		return false;
	}

	std::shared_ptr<E> removeEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		std::shared_ptr<E> e = getEdge(sourceVertex, targetVertex);

		if (e != nullptr) {
			removeEdgeFromTouchingVertices(e);
			EdgeSet.erase(e);
			edgeMap.erase(e);
		}
		return e;
	}

	std::shared_ptr<E> addEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		assertVertexExist(sourceVertex);
		assertVertexExist(targetVertex);
		if (!allowingMultipleEdges && containsEdge(sourceVertex, targetVertex))
			return nullptr;

		if (!allowingLoops && sourceVertex == targetVertex)
			throw std::invalid_argument(LOOPS_NOT_ALLOWED);

		std::shared_ptr<E> e = createEdge(sourceVertex, targetVertex);
		if (containsEdge(e))
			return nullptr;
		edgeMap.insert(std::make_pair(e, IntrusiveEdge<V>(sourceVertex, targetVertex)));
		EdgeSet.insert(e);
		getEdgeContainer(sourceVertex).addOutgoingEdge(e);
		getEdgeContainer(targetVertex).addIncomingEdge(e);
		return e;
	}

	bool containsEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		return getEdge(sourceVertex, targetVertex) != nullptr;
	}

	virtual std::shared_ptr<E>
	createEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		return std::make_shared<E>();
	}

	virtual bool removeVertex(const std::shared_ptr<V> &v) {
		if (containsVertex(v)) {
			removeAllEdges(edgesof(v));
			vertexMapDirected.erase(v);
			VertexSet.erase(v);
			return true;
		}
		return false;
	}

	bool removeAllEdges(const std::vector<std::shared_ptr<E>> &edges) {
		bool modified = false;
		for (const auto &e: edges) {
			modified |= removeEdge(e);
		}

		return modified;
	}

	bool removeAllEdges(const std::list<std::shared_ptr<E>> &edges) {
		bool modified = false;
		for (const auto &e: edges) {
			modified |= removeEdge(e);
		}

		return modified;
	}

	bool removeAllEdges(std::unordered_set<std::shared_ptr<E>> edges) {
		bool modified = false;
		for (const auto &e: edges) {
			modified |= removeEdge(e);
		}

		return modified;
	}

	bool removeAllVertices(const std::vector<std::shared_ptr<V>> &vertices) {
		bool modified = false;

		for (const std::shared_ptr<V> &v: vertices) {
			modified |= removeVertex(v);
		}
		return modified;
	}

	bool removeAllVertices(const std::list<std::shared_ptr<V>> &vertices) {
		bool modified = false;

		for (const std::shared_ptr<V> &v: vertices) {
			modified |= removeVertex(v);
		}
		return modified;
	}

	bool removeAllVertices(std::unordered_set<std::shared_ptr<V>> vertices) {
		bool modified = false;

		for (std::shared_ptr<V> v: vertices) {
			modified |= removeVertex(v);
		}
		return modified;
	}

	bool isRefSink(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");

		for (const std::shared_ptr<E> &e: outgoingEdgesOf(v)) {
			if (e->getIsRef())
				return false;
		}

		for (const std::shared_ptr<E> &e: incomingEdgesOf(v)) {
			if (e->getIsRef())
				return true;
		}

		return VertexSet.size() == 1;
	}

	std::shared_ptr<E> incomingEdgeOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		std::unordered_set<std::shared_ptr<E>> &edgesSet = incomingEdgesOf(v);
		if (edgesSet.size() > 1) {
			throw std::invalid_argument("Cannot get a single incoming edge for a vertex with multiple incoming edges");
		}
		return edgesSet.empty() ? nullptr : *edgesSet.begin();
	}

	std::shared_ptr<E> outgoingEdgeOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		std::unordered_set<std::shared_ptr<E>> &edgesSet = outgoingEdgesOf(v);
		if (edgesSet.size() > 1) {
			throw std::invalid_argument("Cannot get a single incoming edge for a vertex with multiple incoming edges");
		}
		return edgesSet.empty() ? nullptr : *edgesSet.begin();
	}

	std::shared_ptr<V>
	getNextReferenceVertex(const std::shared_ptr<V> &v, bool allowNonRefPaths, std::shared_ptr<E> blacklistedEdge) {
		if (v == nullptr)
			return nullptr;

		for (const std::shared_ptr<E> &edgeToTest: outgoingEdgesOf(v)) {
			if (edgeToTest->getIsRef()) {
				return getEdgeTarget(edgeToTest);
			}
		}

		if (!allowNonRefPaths)
			return nullptr;

		std::vector<std::shared_ptr<E>> edges;
		for (const std::shared_ptr<E> &edgeToTest: outgoingEdgesOf(v)) {
			if (edgeToTest != blacklistedEdge) {
				edges.template emplace_back(edgeToTest);
			}
			if (edges.size() > 2)
				break;
		}
		return edges.size() == 1 ? getEdgeTarget(edges.at(0)) : nullptr;
	}

	std::shared_ptr<V> getPrevReferenceVertex(const std::shared_ptr<V> &v) {
		if (v == nullptr)
			return nullptr;
		std::vector<std::shared_ptr<V>> allVertexs;
		for (const std::shared_ptr<E> &edge: incomingEdgesOf(v)) {
			std::shared_ptr<V> v = getEdgeSource(edge);
			if (isReferenceNode(v))
				allVertexs.template emplace_back(v);
		}
		return allVertexs.size() > 0 ? allVertexs.at(0) : nullptr;
	}

	bool isReferenceNode(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		std::unordered_set<std::shared_ptr<E>> edges = edgesof(v);
		for (const std::shared_ptr<E> &edge: edges) {
			if (edge->getIsRef())
				return true;
		}
		return VertexSet.size() == 1;
	}

	std::shared_ptr<V> getReferenceSourceVertex() {
		for (const std::shared_ptr<V> &vertex: VertexSet) {
			if (isRefSource(vertex))
				return vertex;
		}
		return nullptr;
	}

	std::shared_ptr<V> getReferenceSinkVertex() {
		for (const std::shared_ptr<V> &vertex: VertexSet) {
			if (isRefSink(vertex))
				return vertex;
		}
		return nullptr;
	}


	/**
	 * Remove all vertices in the graph that aren't on a path from the reference source vertex to the reference sink vertex
	 *
	 * More aggressive reference pruning algorithm than removeVerticesNotConnectedToRefRegardlessOfEdgeDirection,
	 * as it requires vertices to not only be connected by a series of directed edges but also prunes away
	 * paths that do not also meet eventually with the reference sink vertex
	 */
	void removePathsNotConnectedToRef(unsigned long size) {
		if (getReferenceSourceVertex() == nullptr || getReferenceSinkVertex() == nullptr) {
			throw std::invalid_argument("Graph must have ref source and sink vertices");
		}
		std::vector<std::shared_ptr<V>> onPathFromRefSource;
		BaseGraphIterator<V, E> sourceIter = BaseGraphIterator<V, E>(this, getReferenceSourceVertex(), false, true);
		//onPathFromRefSource.reserve(size);
		while (sourceIter.hasNext()) {
			onPathFromRefSource.push_back(sourceIter.next());
		}

		std::unordered_set<std::shared_ptr<V>> onPathFromRefSink;
		BaseGraphIterator<V, E> sinkIter = BaseGraphIterator<V, E>(this, getReferenceSinkVertex(), true, false);
		onPathFromRefSink.reserve(size);
		while (sinkIter.hasNext()) {
			onPathFromRefSink.insert(sinkIter.next());
		}

		std::unordered_set<std::shared_ptr<V>> verticesToRemove = getVertexSet();
		for (typename std::vector<std::shared_ptr<V>>::iterator iter = onPathFromRefSource.begin();
		     iter != onPathFromRefSource.end(); iter++) {
			if (onPathFromRefSink.find(*iter) != onPathFromRefSink.end())
				verticesToRemove.erase(*iter);
		}

		removeAllVertices(verticesToRemove);

		if (getSinks().size() > 1)
			throw std::length_error("Should have eliminated all but the reference sink");

		if (getSources().size() > 1)
			throw std::length_error("hould have eliminated all but the reference source");
	}

	std::unordered_set<std::shared_ptr<V>> incomingVerticesOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		std::unordered_set<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<E> &e: incomingEdgesOf(v)) {
			ret.insert(getEdgeSource(e));
		}
		return ret;
	}

	std::unordered_set<std::shared_ptr<V>> outgoingVerticesOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		std::unordered_set<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<E> &e: outgoingEdgesOf(v)) {
			ret.insert(getEdgeTarget(e));
		}
		return ret;
	}

	std::vector<std::shared_ptr<V>> vecIncomingVerticesOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		std::vector<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<E> &e: incomingEdgesOf(v)) {
			ret.push_back(getEdgeSource(e));
		}
		return ret;
	}

	std::vector<std::shared_ptr<V>> vecOutgoingVerticesOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		std::vector<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<E> &e: outgoingEdgesOf(v)) {
			ret.push_back(getEdgeTarget(e));
		}
		return ret;
	}

	std::unordered_set<std::shared_ptr<V>> getSinks() {
		std::unordered_set<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<V> &v: VertexSet) {
			if (isSink(v))
				ret.insert(v);
		}
		return ret;
	}

	std::unordered_set<std::shared_ptr<V>> getSources() {
		std::unordered_set<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<V> &v: VertexSet) {
			if (isSource(v))
				ret.insert(v);
		}
		return ret;
	}

	void cleanNonRefPaths() {
		if (getReferenceSourceVertex() == nullptr || getReferenceSinkVertex() == nullptr)
			return;
		std::unordered_set<std::shared_ptr<E>> edgesToCheck;
		std::shared_ptr<E> e;

		edgesToCheck = incomingEdgesOf(getReferenceSourceVertex());
		while (!edgesToCheck.empty()) {
			e = *(edgesToCheck.begin());
			if (!e->getIsRef()) {
				for (const std::shared_ptr<E> &e: incomingEdgesOf(getEdgeSource(e))) {
					edgesToCheck.insert(e);
				}
				removeEdge(e);
			}
			edgesToCheck.erase(e);
		}

		edgesToCheck = outgoingEdgesOf(getReferenceSinkVertex());
		while (!edgesToCheck.empty()) {
			e = *(edgesToCheck.begin());
			if (!e->getIsRef()) {
				for (const std::shared_ptr<E> &e: outgoingEdgesOf(getEdgeTarget(e))) {
					edgesToCheck.insert(e);
				}
				removeEdge(e);
			}
			edgesToCheck.erase(e);
		}

		Specifics<V, E>::removeSingletonOrphanVertices();
	}


	bool isRefSource(const std::shared_ptr<V> &v) {
		return Specifics<V, E>::isRefSource(v);
	}

	void removeVerticesNotConnectedToRefRegardlessOfEdgeDirection() {
		std::unordered_set<std::shared_ptr<V>> toRemove = VertexSet;
		std::shared_ptr<V> refV = getReferenceSourceVertex();
		if (refV != nullptr) {
			BaseGraphIterator<V, E> iter = BaseGraphIterator<V, E>(this, refV, true, true);
			while (iter.hasNext()) {
				toRemove.erase(iter.next());
			}
		}
		removeAllVertices(toRemove);
	}

	void
	addOrUpdateEdge(const std::shared_ptr<V> &source, const std::shared_ptr<V> &target, const std::shared_ptr<E> &e) {
		if (source.get() == nullptr)
			throw std::invalid_argument("source");
		if (target.get() == nullptr)
			throw std::invalid_argument("target");
		if (e.get() == nullptr)
			throw std::invalid_argument("e");

		std::shared_ptr<E> prev = getEdge(source, target);
		if (prev != nullptr) {
			prev->add(*e);
		} else {
			addEdge(source, target, e);
		}
	}

	std::unordered_set<std::shared_ptr<E>> getEdgeSet() {
		return EdgeSet;
	}

	bool containsAllVertices(std::unordered_set<std::shared_ptr<V>> &vertices) {
		if (vertices.empty())
			throw std::invalid_argument("null vertex");
		for (std::shared_ptr<V> v: vertices) {
			if (!containsVertex(v))
				return false;
		}
		return true;
	}

	virtual void reserveSpace(int size) {
		edgeMap.reserve(2 * size);
		EdgeSet.reserve(2 * size);
		vertexMapDirected.reserve(size);
		VertexSet.reserve(size);
	}

};


#endif //MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H
