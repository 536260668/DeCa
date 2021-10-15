//
// Created by 梦想家xixi on 2021/10/13.
//

#include "ActivityProfileState.h"
#include "Mutect2Utils.h"
#include <iostream>

double ActivityProfileState::isActiveProb() const {return activeProb;}

void ActivityProfileState::setIsActiveProb(const double aP) {this->activeProb = aP;}

Type ActivityProfileState::getResultState() const {return resultState;}

double ActivityProfileState::getResultValue() const {return resultValue;}

ActivityProfileState::ActivityProfileState(SimpleInterval& loc, const double activeProb, const Type resultState, const double resultValue) : loc(loc), activeProb(activeProb), resultState(resultState), resultValue(resultValue){
    Mutect2Utils::validateArg(loc.size() == 1, "Location for an ActivityProfileState must have to size 1 bp.");
    Mutect2Utils::validateArg(resultValue >= 0, "Result value isn't null and its < 0, which is illegal");
}

ActivityProfileState::ActivityProfileState(SimpleInterval& loc, const double activeProb) : loc(loc), activeProb(activeProb), resultState(resultState), resultValue(resultValue){
    Mutect2Utils::validateArg(loc.size() == 1, "Location for an ActivityProfileState must have to size 1 bp.");
    Mutect2Utils::validateArg(resultValue >= 0, "Result value isn't null and its < 0, which is illegal");
}

int ActivityProfileState::getOffset(Locatable *regionStartLoc) {
    Mutect2Utils::validateArg(regionStartLoc != nullptr, "Null object is not allowed here.");
    return loc.getStart() - regionStartLoc->getStart();
}

SimpleInterval ActivityProfileState::getLoc() {
    return loc;
}

std::ostream & operator<<(std::ostream & os, ActivityProfileState&  activityProfileState) {
    std::cout << "loc:  " << activityProfileState.getLoc() << std::endl << "activeProb:" << activityProfileState.isActiveProb() << "   resultValue:"
              << activityProfileState.getResultValue() << "   resultState:" << activityProfileState.getResultState() << std::endl;
    return os;
}

