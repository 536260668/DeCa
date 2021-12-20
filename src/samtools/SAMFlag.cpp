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
