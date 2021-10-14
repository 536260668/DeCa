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

ActivityProfileState::ActivityProfileState(SimpleInterval* loc, const double activeProb, const Type resultState, const double resultValue) {
    Mutect2Utils::validateArg(loc->size() == 1, "Location for an ActivityProfileState must have to size 1 bp.");
    Mutect2Utils::validateArg(resultValue >= 0, "Result value isn't null and its < 0, which is illegal");

    this->loc = loc;
    setIsActiveProb(activeProb);
    this->resultState = resultState;
    this->resultValue = resultValue;
}

ActivityProfileState::ActivityProfileState(SimpleInterval* loc, const double activeProb) {
    new (this)ActivityProfileState(loc, activeProb, NONE, 0);
}

int ActivityProfileState::getOffset(Locatable *regionStartLoc) const {
    Mutect2Utils::validateArg(regionStartLoc != nullptr, "Null object is not allowed here.");
    return loc->getStart() - regionStartLoc->getStart();
}

SimpleInterval* ActivityProfileState::getLoc() const {
    return loc;
}

std::ostream & operator<<(std::ostream & os, const ActivityProfileState * activityProfileState) {
    std::cout << "loc:  " << activityProfileState->loc << std::endl << "activeProb:" << activityProfileState->activeProb << "   resultValue:"
        << activityProfileState->resultValue << "   resultState:" << activityProfileState->resultState << std::endl;
    return os;
}

