//
// Created by 梦想家xixi on 2021/10/21.
//

#ifndef MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H
#define MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H

#include "parallel_hashmap/phmap.h"
#include <memory>

template<class EE>
class DirectedEdgeContainer {
public:
    phmap::flat_hash_set<std::shared_ptr<EE>> incoming;

    phmap::flat_hash_set<std::shared_ptr<EE>> outgoing;

	void addIncomingEdge(const std::shared_ptr<EE> &e) { incoming.insert(e); }

	void addOutgoingEdge(const std::shared_ptr<EE> &e) { outgoing.insert(e); }

    phmap::flat_hash_set<std::shared_ptr<EE>> &getUnmodifiableIncomingEdges() {
		return incoming;
	}

    phmap::flat_hash_set<std::shared_ptr<EE>> &getUnmodifiableOutgoingEdges() {
		return outgoing;
	}

	void removeIncomingEdge(const std::shared_ptr<EE> &e) { incoming.erase(e); }

	void removeOutgoingEdge(const std::shared_ptr<EE> &e) { outgoing.erase(e); }
};


#endif //MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H
