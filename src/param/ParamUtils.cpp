//
// Created by 梦想家xixi on 2021/12/3.
//

#include "ParamUtils.h"
#include <utility>


int ParamUtils::isPositiveOrZero(int val, std::string message) {
    Mutect2Utils::validateArg(val >= 0, std::move(message));
    return val;
}
