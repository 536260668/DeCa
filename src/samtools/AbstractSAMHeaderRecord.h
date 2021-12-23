//
// Created by 梦想家xixi on 2021/12/22.
//

#ifndef MUTECT2CPP_MASTER_ABSTRACTSAMHEADERRECORD_H
#define MUTECT2CPP_MASTER_ABSTRACTSAMHEADERRECORD_H
#include <map>
#include <string>


class AbstractSAMHeaderRecord {
private:
    std::map<std::string, std::string> mAttributes;
public:
    AbstractSAMHeaderRecord() = default;
    virtual ~AbstractSAMHeaderRecord() = 0;
    std::string getAttribute(std::string & key);
    virtual void setAttribute(std::string & key, std::string & value);
    std::map<std::string, std::string>& getAttributes();
};


#endif //MUTECT2CPP_MASTER_ABSTRACTSAMHEADERRECORD_H
