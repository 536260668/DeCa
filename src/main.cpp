//
// Created by 梦想家xixi on 2021/10/11.
//

#include <iostream>
#include "ActivityProfileState.h"

int main(int argc, char **argv)
{
    SimpleInterval* simpleInterval = new SimpleInterval("chr1", 100, 100);
    SimpleInterval* simpleInterval2 = new SimpleInterval("chr1", 130, 170);
    ActivityProfileState* activityProfileState = new ActivityProfileState(simpleInterval, 0.11);
    std::cout << activityProfileState << std::endl;
    std::cout << activityProfileState->getOffset(simpleInterval2) << std::endl;
}
