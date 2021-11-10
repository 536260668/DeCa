//
// Created by 梦想家xixi on 2021/11/8.
//

#include "StringUtils.h"

void StringUtils::toUpperCase(uint8_t* & bytes, int length) {
    for(int i = 0; i < length; ++i) {
        if (bytes[i] >= 97 && bytes[i] <= 122) {
            bytes[i] += -32;
        }
    }
}
