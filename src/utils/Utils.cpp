//
// Created by lhh on 4/26/22.
//

#include <cstring>
#include "Utils.h"

bool Utils::equalRange(uint8_t *left, int leftOffset, uint8_t *right, int rightOffset, int length)
{
    for(int i=0; i<length; i++)
    {
        if(left[leftOffset + i] != right[rightOffset + i])
            return false;
    }
    return true;
}

std::shared_ptr<char[]> Utils::dupBytes(char b, int nCopies)
{
    std::shared_ptr<char[]> bytes(new char[nCopies]);
    memset(bytes.get(), (int)b, nCopies);
    return bytes;
}