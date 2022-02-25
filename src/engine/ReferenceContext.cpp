//
// Created by 梦想家xixi on 2022/1/4.
//

#include "ReferenceContext.h"

#include <utility>

ReferenceContext::ReferenceContext(char *dataSource, std::shared_ptr<SimpleInterval>  interval) : dataSource(dataSource), interval(std::move(interval)){

}

uint8_t ReferenceContext::getBase() {
    return dataSource[interval->getStart()];
}

const std::shared_ptr<SimpleInterval> & ReferenceContext::getInterval() {
    return interval;
}
