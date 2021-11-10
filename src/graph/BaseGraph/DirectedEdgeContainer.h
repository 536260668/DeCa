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
    ArraySet<EE*> unmodifiableIncoming;
    ArraySet<EE*> unmodifiableOutgoing;

public:
    ArraySet<EE*> incoming;
    ArraySet<EE*> outgoing;

    void addIncomingEdge(EE* e) {incoming.insert(e);}

    void addOutgoingEdge(EE* e) {outgoing.insert(e);}

    ArraySet<EE*>  getUnmodifiableIncomingEdges() {
        unmodifiableIncoming = incoming;
        return unmodifiableIncoming;
    }

    ArraySet<EE*> getUnmodifiableOutgoingEdges() {
        unmodifiableOutgoing = outgoing;
        return unmodifiableOutgoing;
    }

    void removeIncomingEdge(EE* e) {
        incoming.erase(e);
    }

    void removeOutgoingEdge(EE* e) {
        outgoing.erase(e);
    }
};



#endif //MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H
