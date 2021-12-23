//
// Created by 梦想家xixi on 2021/12/22.
//

#include "AbstractSAMHeaderRecord.h"

AbstractSAMHeaderRecord::~AbstractSAMHeaderRecord() = default;

std::string AbstractSAMHeaderRecord::getAttribute(std::string &key) {
    return mAttributes.at(key);
}

void AbstractSAMHeaderRecord::setAttribute(std::string &key, std::string &value) {
    if(value.empty()) {
        mAttributes.erase(key);
    } else {
        mAttributes.insert(std::pair<std::string, std::string>(key, value));
    }
}

std::map<std::string, std::string>& AbstractSAMHeaderRecord::getAttributes() {
    return mAttributes;
}



