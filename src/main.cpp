//
// Created by 梦想家xixi on 2021/10/11.
//

#include <iostream>
#include "ActivityProfileState.h"

int main(int argc, char **argv)
{
    SimpleInterval simpleInterval1("chr1", 100, 200);
    SimpleInterval simpleInterval2("chr1", 150, 250);
    SimpleInterval* ptr = simpleInterval1.expandWithinContig(90,20);
    std::cout << *ptr << std::endl;
    delete ptr;
}
