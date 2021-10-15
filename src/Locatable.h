//
// Created by 梦想家xixi on 2021/10/12.
//

#ifndef MUTECT2CPP_MASTER_LOCATABLE_H
#define MUTECT2CPP_MASTER_LOCATABLE_H

#include "Mutect2Utils.h"
#include <string>

class Locatable {
public:
    virtual std::string getContig() const = 0;
    virtual int getStart() const = 0;
    virtual int getEnd() const = 0;

    virtual int getLengthOnReference();
    virtual bool contains(Locatable* other);
    virtual bool contigsMatch(Locatable* other);
    virtual bool withinDistanceOf(Locatable* other, int distance);
    virtual bool overlaps(Locatable* other);
};


#endif //MUTECT2CPP_MASTER_LOCATABLE_H
