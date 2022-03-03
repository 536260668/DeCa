//
// Created by 梦想家xixi on 2021/10/21.
//

#ifndef MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H
#define MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H

#include <unordered_set>

template<class VV, class EE>
class DirectedEdgeContainer {
private:
    static const long serialVersionUID = 7494242245729767106L;
//    std::set<std::shared_ptr<EE>> unmodifiableIncoming;
//    std::set<std::shared_ptr<EE>> unmodifiableOutgoing;

public:
    std::unordered_set<std::shared_ptr<EE>> incoming;
    std::unordered_set<std::shared_ptr<EE>> outgoing;

    void addIncomingEdge(const std::shared_ptr<EE> & e) {incoming.insert(e);}

    void addOutgoingEdge(const std::shared_ptr<EE> & e) {outgoing.insert(e);}

    std::unordered_set<std::shared_ptr<EE>>& getUnmodifiableIncomingEdges() {
        //unmodifiableIncoming = incoming;
        return incoming;
    }

    std::unordered_set<std::shared_ptr<EE>>& getUnmodifiableOutgoingEdges() {
        //unmodifiableOutgoing = outgoing;
        return outgoing;
    }

    void removeIncomingEdge(const std::shared_ptr<EE> & e) {
        incoming.erase(e);
    }

    void removeOutgoingEdge(const std::shared_ptr<EE> & e) {
        outgoing.erase(e);
    }
};



#endif //MUTECT2CPP_MASTER_DIRECTEDEDGECONTAINER_H
