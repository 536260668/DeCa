//
// Created by 梦想家xixi on 2021/12/22.
//

#include "SAMFileHeader.h"

const std::string SAMFileHeader::VERSION_TAG = "VN";
const std::string SAMFileHeader::SORT_ORDER_TAG = "SO";
const std::string SAMFileHeader::GROUP_ORDER_TAG = "GO";
const std::string SAMFileHeader::CURRENT_VERSION = "1.6";

const std::set<std::string> SAMFileHeader::ACCEPTABLE_VERSIONS{"1.0", "1.3", "1.4", "1.5", "1.6"};
const std::set<std::string> SAMFileHeader::STANDARD_TAGS{"VN", "SO", "GO"};

SAMFileHeader::SAMFileHeader() {
}

void SAMFileHeader::init() {
    setAttribute((std::string &) "VN", (std::string &) "1.6");
}

void SAMFileHeader::setAttribute(std::string &key, std::string &value) {
    std::string tempVal = value;
    if(key == "SO") {
        sortOrder = SAMFileHeader_SortOrder_null;
        if(value == "SAMFileHeader_unsorted" ||value == "SAMFileHeader_queryname" ||value == "SAMFileHeader_coordinate" ||value == "SAMFileHeader_duplicate" ||value == "SAMFileHeader_SortOrder_null") {
            tempVal = value;
        } else {
            tempVal = "SAMFileHeader_unknown";
        }
    } else if (key == "GO") {
        groupOrder = SAMFileHeader_GroupOrder_null;
    }
    AbstractSAMHeaderRecord::setAttribute(key, tempVal);
}

int SAMFileHeader::getSequenceIndex(std::string &basicString) {
    return mSequenceDictionary.getSequenceIndex(basicString);
}
