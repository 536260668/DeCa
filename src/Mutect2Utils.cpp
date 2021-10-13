//
// Created by 梦想家xixi on 2021/10/12.
//
#include <iostream>
#include <stdexcept>

#include "Mutect2Utils.h"

std::string Mutect2Utils::replaceWith(std::string& str1, const std::string& str2, const std::string& str3){
    int pos;
    pos = str1.find(str2);
    while(pos != -1){
        str1.replace(pos,str1.length(),str3);
        pos = str1.find(str2);
    }
    std::cout << str1 << std::endl;
    return str1;
}

bool Mutect2Utils::overlaps(int start, int end, int start2, int end2) {
    return start2 >= start && start2 <= end || end2 >= start && end2 <= end || encloses(start2, end2, start, end);
}

bool Mutect2Utils::encloses(int outerStart, int outerEnd, int innerStart, int innerEnd) {
    return innerStart >= outerStart && innerEnd <= outerEnd;
}

void Mutect2Utils::validateArg(bool condition, std::string msg) {
    if(!condition){
        throw std::invalid_argument(msg);
    }
}
