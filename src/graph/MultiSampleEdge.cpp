//
// Created by 梦想家xixi on 2021/10/20.
//

#include "MultiSampleEdge.h"
#include "Mutect2Utils.h"

MultiSampleEdge::MultiSampleEdge(bool isRef, int multiplicity, int singleSampleCapacity) : BaseEdge(isRef, multiplicity), singleSampleCapacity(singleSampleCapacity){
    Mutect2Utils::validateArg(singleSampleCapacity > 0, "singleSampleCapacity must be > 0");
    singleSampleMultiplicities.push(singleSampleCapacity);
    singleSampleMultiplicities.push(multiplicity);
}

void MultiSampleEdge::flushSingleSampleMultiplicity() {
    singleSampleMultiplicities.push(currentSingleSampleMultiplicity);
    if(singleSampleMultiplicities.size() == singleSampleCapacity + 1) {
        singleSampleMultiplicities.pop();
    } else if (singleSampleMultiplicities.size() > singleSampleCapacity + 1) {
        throw std::invalid_argument("Somehow the per sample multiplicity list has grown too big");
    }
    currentSingleSampleMultiplicity = 0;
}

void MultiSampleEdge::incMultiplicity(int incr) {
    BaseEdge::incMultiplicity(incr);
    currentSingleSampleMultiplicity += incr;
}

int MultiSampleEdge::getPruningMultiplicity() const{
    return singleSampleMultiplicities.top();
}
