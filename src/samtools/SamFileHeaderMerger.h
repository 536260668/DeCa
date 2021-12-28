//
// Created by 梦想家xixi on 2021/12/27.
//

#ifndef MUTECT2CPP_MASTER_SAMFILEHEADERMERGER_H
#define MUTECT2CPP_MASTER_SAMFILEHEADERMERGER_H

#include "SAMFileHeader.h"
#include "HeaderRecordAndFileHeader.h"

class SamFileHeaderMerger {
private:
    int recordCounter;
    static char INT_TO_BASE36[36];
    bool hasProgramGroupCollisions;
    bool hasReadGroupCollisions;
    std::vector<SAMProgramRecord*> pgToDelete;
    std::vector<SAMReadGroupRecord*> rgToDelete;
    std::map<SAMFileHeader*, std::map<std::string, std::string>> samProgramGroupIdTranslation;

    bool mergeHeaderRecords(std::vector<HeaderRecordAndFileHeader> &headerRecords, HeaderRecordFactory* headerRecordFactory, std::set<std::string> & idsThatAreAlreadyTaken, std::map<SAMFileHeader*, std::map<std::string, std::string>> & idTranslationTable,
                            std::vector<AbstractSAMHeaderRecord*> & result);

    std::vector<SAMProgramRecord> mergeProgramGroups(std::vector<SAMFileHeader*> headers);

    std::vector<HeaderRecordAndFileHeader> translateIds(const std::vector<HeaderRecordAndFileHeader>& programGroups, std::map<SAMFileHeader*, std::map<std::string, std::string>>, std::vector<SAMProgramRecord*> &todelete, bool translatePpIds);

    static bool compareById (AbstractSAMHeaderRecord& a, AbstractSAMHeaderRecord& b) {return a.getId() < b.getId();}

public:
    std::string positiveFourDigitBase36Str(int leftOver);

    SamFileHeaderMerger();
};


#endif //MUTECT2CPP_MASTER_SAMFILEHEADERMERGER_H
