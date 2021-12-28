//
// Created by 梦想家xixi on 2021/12/27.
//

#ifndef MUTECT2CPP_MASTER_HEADERRECORDFACTORY_H
#define MUTECT2CPP_MASTER_HEADERRECORDFACTORY_H

#include "AbstractSAMHeaderRecord.h"

class HeaderRecordFactory {
public:
    virtual AbstractSAMHeaderRecord* createRecord(std::string & newId, AbstractSAMHeaderRecord* record);
};


#endif //MUTECT2CPP_MASTER_HEADERRECORDFACTORY_H
