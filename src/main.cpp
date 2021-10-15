//
// Created by 梦想家xixi on 2021/10/11.
//

#include <iostream>
#include "ActivityProfileState.h"
#include "AssemblyRegion.h"

int main(int argc, char **argv)
{
    SimpleInterval s1("chr1", 1, 1);
    SimpleInterval s2("chr1", 2, 2);
    SimpleInterval s3("chr1", 3, 3);
    SimpleInterval s4("chr1", 1, 3);
    ActivityProfileState a1(s1, 0.1);
    ActivityProfileState a2(s2, 0.2);
    ActivityProfileState a3(s3, 0.3);
    std::vector<ActivityProfileState> ac;
    ac.push_back(a1);
    ac.push_back(a2);
    ac.push_back(a3);
    AssemblyRegion as(s4, ac, true, 1);
    std::cout << as << std::endl;
}
