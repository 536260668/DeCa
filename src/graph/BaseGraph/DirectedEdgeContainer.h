//
// Created by 梦想家xixi on 2021/10/21.
//

#ifndef MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H
#define MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H

#include "../set/ArraySet.h"

template<class VV, class EE>
class DirectedEdgeContainer {
private:
    static const long serialVersionUID = 7494242245729767106L;
    ArraySet<std::shared_ptr<EE>> unmodifiableIncoming;
    ArraySet<std::shared_ptr<EE>> unmodifiableOutgoing;

public:
    ArraySet<std::shared_ptr<EE>> incoming;
    ArraySet<std::shared_ptr<EE>> outgoing;

    void addIncomingEdge(std::shared_ptr<EE> e) {incoming.insert(e);}

    void addOutgoingEdge(std::shared_ptr<EE> e) {outgoing.insert(e);}

    ArraySet<std::shared_ptr<EE>>  getUnmodifiableIncomingEdges() {
        unmodifiableIncoming = incoming;
        return unmodifiableIncoming;
    }

    ArraySet<std::shared_ptr<EE>> getUnmodifiableOutgoingEdges() {
        unmodifiableOutgoing = outgoing;
        return unmodifiableOutgoing;
    }

    void removeIncomingEdge(std::shared_ptr<EE> e) {
        incoming.erase(e);
    }

    void removeOutgoingEdge(std::shared_ptr<EE> e) {
        outgoing.erase(e);
    }
};



#endif //MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H
