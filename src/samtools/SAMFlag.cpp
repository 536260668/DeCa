//
// Created by 梦想家xixi on 2021/12/20.
//

#include "SAMFlag.h"

SAMFlag::SAMFlag(int flag, std::string& description) : flag(flag), description(description){

}

SAMFlag &SAMFlag::READ_PAIRED() {
    static SAMFlag ret(1, (std::string &) "Template having multiple segments in sequencing");
    return ret;
}

SAMFlag &SAMFlag::READ_UNMAPPED() {
    static SAMFlag ret(4, (std::string &) "Segment unmapped");
    return ret;
}

SAMFlag &SAMFlag::MATE_UNMAPPED() {
    static SAMFlag ret(8, (std::string &) "Next segment in the template unmapped");
    return ret;
}

SAMFlag &SAMFlag::PROPER_PAIR() {
    static SAMFlag ret(2, (std::string &) "Each segment properly aligned according to the aligner");
    return ret;
}

SAMFlag &SAMFlag::READ_REVERSE_STRAND() {
    static SAMFlag ret(16, (std::string &) "SEQ being reverse complemented");
    return ret;
}

SAMFlag &SAMFlag::MATE_REVERSE_STRAND() {
    static SAMFlag ret(32, (std::string &) "SEQ of the next segment in the template being reverse complemented");
    return ret;
}

SAMFlag &SAMFlag::FIRST_OF_PAIR() {
    static SAMFlag ret(64, (std::string &) "The first segment in the template");
    return ret;
}

SAMFlag &SAMFlag::SECOND_OF_PAIR() {
    static SAMFlag ret(128, (std::string &) "The last segment in the template");
    return ret;
}

SAMFlag &SAMFlag::SECONDARY_ALIGNMENT() {
    static SAMFlag ret(256, (std::string &) "Secondary alignment");
    return ret;
}

SAMFlag &SAMFlag::READ_FAILS_VENDOR_QUALITY_CHECK() {
    static SAMFlag ret(512, (std::string &) "Not passing quality controls");
    return ret;
}

SAMFlag &SAMFlag::DUPLICATE_READ() {
    static SAMFlag ret(1024, (std::string &) "PCR or optical duplicate");
    return ret;
}

SAMFlag &SAMFlag::SUPPLEMENTARY_ALIGNMENT() {
    static SAMFlag ret(2048, (std::string &) "Supplementary alignment");
    return ret;
}
