//
// Created by 梦想家xixi on 2021/12/22.
//

#include "AbstractSAMHeaderRecord.h"

AbstractSAMHeaderRecord::~AbstractSAMHeaderRecord() = default;

std::string AbstractSAMHeaderRecord::getAttribute(std::string key) {
    if(mAttributes.find(key) != mAttributes.end())
        return mAttributes.at(key);
    else
        return "";
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

std::string & AbstractSAMHeaderRecord::getId() {
    throw std::invalid_argument("Method not implemented");
}



