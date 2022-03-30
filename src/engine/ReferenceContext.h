//
// Created by 梦想家xixi on 2022/1/4.
//

#ifndef MUTECT2CPP_MASTER_REFERENCECONTEXT_H
#define MUTECT2CPP_MASTER_REFERENCECONTEXT_H


#include "SimpleInterval.h"

class ReferenceContext {
private:
    char refBase;
    std::shared_ptr<SimpleInterval> interval;

public:
    ReferenceContext(std::shared_ptr<SimpleInterval>  interval, char refBase);
    uint8_t getBase();
    const std::shared_ptr<SimpleInterval> & getInterval();
};


#endif //MUTECT2CPP_MASTER_REFERENCECONTEXT_H
