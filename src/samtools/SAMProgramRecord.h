//
// Created by 梦想家xixi on 2021/12/22.
//

#ifndef MUTECT2CPP_MASTER_SAMPROGRAMRECORD_H
#define MUTECT2CPP_MASTER_SAMPROGRAMRECORD_H

#include "AbstractSAMHeaderRecord.h"

class SAMProgramRecord : public AbstractSAMHeaderRecord{
private:
    std::string mProgramGroupId;

public:
    static const std::string PROGRAM_GROUP_ID_TAG;
    static const std::string PROGRAM_NAME_TAG;
    static const std::string PROGRAM_VERSION_TAG;
    static const std::string COMMAND_LINE_TAG;
    static const std::string PREVIOUS_PROGRAM_GROUP_ID_TAG;
    SAMProgramRecord(std::string & programGroupId);
    SAMProgramRecord(std::string & id, SAMProgramRecord & srcProgramRecord);
};


#endif //MUTECT2CPP_MASTER_SAMPROGRAMRECORD_H
