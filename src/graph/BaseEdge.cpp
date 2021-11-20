//
// Created by 梦想家xixi on 2021/10/20.
//

#include "BaseEdge.h"
#include "Mutect2Utils.h"

BaseEdge::BaseEdge(bool isRef, int multiplicity) : isRef(isRef), multiplicity(multiplicity){
    Mutect2Utils::validateArg(multiplicity >= 0, "multiplicity must be >= 0");
}

void BaseEdge::incMultiplicity(const int incr) {
    Mutect2Utils::validateArg(incr >= 0, "incr must be >= 0");
    multiplicity += incr;
}

void BaseEdge::setMultiplicity(int value) {
    Mutect2Utils::validateArg(value >= 0, "multiplicity must be >= 0");
    multiplicity = value;
}

void BaseEdge::setIsRef(const bool ref) {
    this->isRef = ref;
}

BaseEdge BaseEdge::add(BaseEdge &edge) {
    multiplicity += edge.getMultiplicity();
    isRef = isRef || edge.getIsRef();
    return *this;
}

BaseEdge *BaseEdge::makeOREdge(std::vector<BaseEdge *> edges, int multiplicity) {
    Mutect2Utils::validateArg(!edges.empty(), "have no edge");
    bool anyRef = false;
    for(BaseEdge* edge : edges) {
        if(edge->getIsRef()) {
            anyRef = true;
            break;
        }
    }
    return new BaseEdge(anyRef, multiplicity);
}

