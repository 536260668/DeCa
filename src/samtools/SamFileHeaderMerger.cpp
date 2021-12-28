//
// Created by 梦想家xixi on 2021/12/27.
//

#include "SamFileHeaderMerger.h"
#include <vector>

char SamFileHeaderMerger::INT_TO_BASE36[36]{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};



bool SamFileHeaderMerger::mergeHeaderRecords(std::vector<HeaderRecordAndFileHeader> &headerRecords, HeaderRecordFactory* headerRecordFactory,
                                             std::set<std::string> &idsThatAreAlreadyTaken,
                                             std::map<SAMFileHeader*, std::map<std::string, std::string>> &idTranslationTable,
                                             std::vector<AbstractSAMHeaderRecord*> &result) {
    std::map<std::string, std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>>> idToRecord;
    for(HeaderRecordAndFileHeader& pair : headerRecords) {
        AbstractSAMHeaderRecord* record = pair.getHeaderRecord();
        SAMFileHeader* header = pair.getFileHeader();
        std::string recordId = record->getId();
        if(idToRecord.find(recordId) == idToRecord.end()) {
            idToRecord.insert(std::pair<std::string, std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>>>(recordId, std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>>()));
        }
        std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>>& recordsWithSameId = idToRecord.at(recordId);

        if(recordsWithSameId.find(record) == recordsWithSameId.end()) {
            recordsWithSameId.insert(std::pair<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>>(record, std::vector<SAMFileHeader*>()));
        }
        std::vector<SAMFileHeader*>& fileHeaders = recordsWithSameId.at(record);
        fileHeaders.emplace_back(header);
    }
    bool hasCollisions = false;
    std::map<std::string, std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>>>::iterator miter;
    for(miter = idToRecord.begin(); miter != idToRecord.end(); miter++) {
        std::string recordId = miter->first;
        std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>> & recordsWithSameId = miter->second;
        std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>>::iterator iter;
        for(iter = recordsWithSameId.begin(); iter != recordsWithSameId.end(); iter++) {
            AbstractSAMHeaderRecord* record = iter->first;
            std::vector<SAMFileHeader*> & fileHeaders = iter->second;
            std::string newId;
            if(idsThatAreAlreadyTaken.find(recordId) == idsThatAreAlreadyTaken.end()) {
                newId = recordId;
                idsThatAreAlreadyTaken.insert(recordId);
                ++recordCounter;
            } else {
                hasCollisions = true;
                while(idsThatAreAlreadyTaken.find(newId = recordId + "." + positiveFourDigitBase36Str(recordCounter++)) != idsThatAreAlreadyTaken.end());
                idsThatAreAlreadyTaken.insert(newId);
            }
            for(SAMFileHeader* fileHeader : fileHeaders) {
                if(idTranslationTable.find(fileHeader) == idTranslationTable.end()) {
                    idTranslationTable.insert(std::pair<SAMFileHeader*, std::map<std::string, std::string>>(fileHeader, std::map<std::string, std::string>()));
                }
                std::map<std::string, std::string> & readerTranslationTable = idTranslationTable.at(fileHeader);
                readerTranslationTable.insert(std::pair<std::string, std::string>(newId, recordId));
            }
            result.emplace_back(headerRecordFactory->createRecord(newId, record));
        }
    }
    return hasCollisions;
}

std::string SamFileHeaderMerger::positiveFourDigitBase36Str(int leftOver) {
    if(leftOver == 0) {
        return "0";
    } else {
        std::string ret;
        while (leftOver > 0) {
            int valueIndex = leftOver % 36;
            ret += INT_TO_BASE36[valueIndex];
            leftOver /= 36;
        }
        std::reverse(ret.begin(), ret.end());
        return ret;
    }
}

std::vector<SAMProgramRecord> SamFileHeaderMerger::mergeProgramGroups(std::vector<SAMFileHeader *> headers) {
    std::vector<SAMProgramRecord> overallResult;
    std::set<std::string> idsThatAreAlreadyTaken;
    std::vector<HeaderRecordAndFileHeader> programGroupsLeftToProcess;
    for(SAMFileHeader* header : headers) {
        for(SAMProgramRecord& programRecord : header->getProgramRecords()) {
            if(idsThatAreAlreadyTaken.find(programRecord.getId()) != idsThatAreAlreadyTaken.end()) {
                throw std::invalid_argument("Input file contains more than one PG with the same id");
            } else {
                idsThatAreAlreadyTaken.insert(programRecord.getId());
            }
            programGroupsLeftToProcess.emplace_back(HeaderRecordAndFileHeader(header, &programRecord));
        }
        idsThatAreAlreadyTaken.clear();
    }
    recordCounter = 0;
    std::vector<HeaderRecordAndFileHeader> currentProgramGroups;
    for(std::vector<HeaderRecordAndFileHeader>::iterator programGroupsLeftToProcessIterator = programGroupsLeftToProcess.begin();
    programGroupsLeftToProcessIterator != programGroupsLeftToProcess.end(); ) {
        HeaderRecordAndFileHeader pair = *programGroupsLeftToProcessIterator;
        if(pair.getHeaderRecord()->getAttribute(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG).empty()) {
            programGroupsLeftToProcess.erase(programGroupsLeftToProcessIterator);
            currentProgramGroups.emplace_back(pair);
        } else {
            programGroupsLeftToProcessIterator++;
        }
    }
    while(!currentProgramGroups.empty()) {
        std::vector<AbstractSAMHeaderRecord*> currentResult;
        std::string tmp = "tmp";
        HeaderRecordFactory* factory = new SAMProgramRecord(tmp);
        hasProgramGroupCollisions |= mergeHeaderRecords(currentProgramGroups, factory, idsThatAreAlreadyTaken, samProgramGroupIdTranslation, currentResult);
        for(AbstractSAMHeaderRecord* samHeaderRecord : currentResult) {
            overallResult.emplace_back(*(SAMProgramRecord*)samHeaderRecord);
        }
        currentProgramGroups = translateIds(currentProgramGroups, samProgramGroupIdTranslation, pgToDelete, false);
        programGroupsLeftToProcess = translateIds(programGroupsLeftToProcess, samProgramGroupIdTranslation, pgToDelete, true);

        std::vector<HeaderRecordAndFileHeader> programGroupsToProcessNext;
        for(std::vector<HeaderRecordAndFileHeader>::iterator programGroupsLeftToProcessIterator = programGroupsLeftToProcess.begin();
        programGroupsLeftToProcessIterator != programGroupsLeftToProcess.end();) {
            std::vector<HeaderRecordAndFileHeader>::iterator tmp = programGroupsLeftToProcessIterator;
            std::string ppIdOfRecordLeftToProcess = programGroupsLeftToProcessIterator->getHeaderRecord()->getAttribute(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG);
            for(HeaderRecordAndFileHeader justProcessedPair : currentProgramGroups) {
                std::string idJustProcessed = justProcessedPair.getHeaderRecord()->getId();
                if(programGroupsLeftToProcessIterator->getFileHeader() == justProcessedPair.getFileHeader() && ppIdOfRecordLeftToProcess == idJustProcessed) {
                    programGroupsLeftToProcess.erase(programGroupsLeftToProcessIterator);
                    programGroupsToProcessNext.emplace_back(*programGroupsLeftToProcessIterator);
                    break;
                }
            }
            if(tmp == programGroupsLeftToProcessIterator) {
                programGroupsLeftToProcessIterator++;
            }
        }
        currentProgramGroups = programGroupsToProcessNext;
    }
    if(!programGroupsLeftToProcess.empty()) {
        throw std::invalid_argument("program groups weren't processed. Do their PP ids point to existing PGs?");
    }
    std::sort(overallResult.begin(), overallResult.end(), compareById);
    return overallResult;
}

std::vector<HeaderRecordAndFileHeader> SamFileHeaderMerger::translateIds(const std::vector<HeaderRecordAndFileHeader>& programGroups,
                                  std::map<SAMFileHeader *, std::map<std::string, std::string>> idTranslationTable, std::vector<SAMProgramRecord*> &todelete, bool translatePpIds) {
    std::vector<HeaderRecordAndFileHeader> result;
    for(HeaderRecordAndFileHeader pair : programGroups) {
        SAMProgramRecord& record = *(SAMProgramRecord*)pair.getHeaderRecord();
        std::string id = record.getId();
        std::string ppId = record.getAttribute(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG);
        SAMFileHeader* header = pair.getFileHeader();
        std::map<std::string, std::string> translations = idTranslationTable.find(header) != idTranslationTable.end() ? idTranslationTable.at(header) : std::map<std::string, std::string>();
        SAMProgramRecord* translatedRecord = nullptr;
        if(!translations.empty()) {
            std::string translatedId = translations.find(id) != translations.end() ? translations.at(id) : "";
            std::string translatedPpId = translatePpIds ? (translations.find(ppId) != translations.end() ? translations.at(ppId) : "") : "";
            bool needToTranslateId = !translatedId.empty() && translatedId != id;
            bool needToTranslatePpId = !translatedPpId.empty() && translatedPpId != ppId;
            if(needToTranslateId && needToTranslatePpId) {
                translatedRecord = new SAMProgramRecord(translatedId, record);
                translatedRecord->setAttribute(
                        const_cast<std::string &>(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG), translatedPpId);
            } else if (needToTranslateId) {
                translatedRecord = new SAMProgramRecord(translatedId, record);
            } else if (needToTranslatePpId) {
                translatedRecord = new SAMProgramRecord(id, record);
                translatedRecord->setAttribute(
                        const_cast<std::string &>(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG), translatedPpId);
            }
        }
        if(translatedRecord != nullptr) {
            result.emplace_back(HeaderRecordAndFileHeader(header, translatedRecord));
            todelete.emplace_back(translatedRecord);
        } else {
            result.emplace_back(pair);
        }
    }
    return result;
}

SamFileHeaderMerger::SamFileHeaderMerger() {
    hasProgramGroupCollisions = false;
    hasReadGroupCollisions = false;
}
