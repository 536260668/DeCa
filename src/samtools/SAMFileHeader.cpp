//
// Created by 梦想家xixi on 2021/12/22.
//

#include "SAMFileHeader.h"

const std::string SAMFileHeader::READ_GROUP_ID_TAG = "ID";
const std::string SAMFileHeader::SEQUENCING_CENTER_TAG = "CN";
const std::string SAMFileHeader::DESCRIPTION_TAG = "DS";
const std::string SAMFileHeader::DATE_RUN_PRODUCED_TAG = "DT";
const std::string SAMFileHeader::FLOW_ORDER_TAG = "FO";
const std::string SAMFileHeader::KEY_SEQUENCE_TAG = "KS";
const std::string SAMFileHeader::LIBRARY_TAG = "LB";
const std::string SAMFileHeader::PROGRAM_GROUP_TAG = "PG";
const std::string SAMFileHeader::PREDICTED_MEDIAN_INSERT_SIZE_TAG = "PI";
const std::string SAMFileHeader::PLATFORM_TAG = "PL";
const std::string SAMFileHeader::PLATFORM_MODEL_TAG = "PM";
const std::string SAMFileHeader::PLATFORM_UNIT_TAG = "PU";
const std::string SAMFileHeader::READ_GROUP_SAMPLE_TAG = "SM";
const std::string SAMFileHeader::BARCODE_TAG = "BC";
const std::set<std::string> SAMFileHeader::STANDARD_TAGS{"ID", "CN", "DS", "DT", "FO", "KS", "LB", "PG", "PI", "PL", "PM", "PU", "SM", "BC"};