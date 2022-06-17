//
// Created by 梦想家xixi on 2022/1/4.
//

#include "ReferenceContext.h"
#include <utility>

ReferenceContext::ReferenceContext(std::shared_ptr<SimpleInterval>  interval, char refBase) : interval(std::move(interval)), refBase(refBase){

}

ReferenceContext::ReferenceContext(const ReferenceContext &other) : interval(std::move(other.interval)), refBase(other.refBase) {
}

uint8_t ReferenceContext::getBase() {
    return refBase;
}

const std::shared_ptr<SimpleInterval> & ReferenceContext::getInterval() {
    return interval;
}
