//
// Created by 梦想家xixi on 2022/1/4.
//

#include "ReferenceContext.h"

ReferenceContext::ReferenceContext(char *dataSource, SimpleInterval &interval) : dataSource(dataSource), interval(interval){

}

uint8_t ReferenceContext::getBase() {
    return dataSource[interval.getStart()];
}

SimpleInterval &ReferenceContext::getInterval() {
    return interval;
}
