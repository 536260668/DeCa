//
// Created by 梦想家xixi on 2021/10/11.
//

#include <iostream>
#include "SimpleInterval.h"

int main(int argc, char **argv)
{
    SimpleInterval* interval = new SimpleInterval("chr1", 1, 100);
    SimpleInterval* interval2 = new SimpleInterval("chr1", 98, 120);
    std::cout << interval->overlaps(interval2) << std::endl;
    interval= interval->intersect(interval2);
    interval->printfInterval();
    delete interval;
    return 0;
}
