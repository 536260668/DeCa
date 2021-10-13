//
// Created by 梦想家xixi on 2021/10/11.
//

#include "SimpleInterval.h"


SimpleInterval::SimpleInterval(std::string str) {
    Mutect2Utils::validateArg(!str.empty(), "Null object is not allowed here.");

    std::string contig;
    int start;
    int end;

    const int colonIndex = str.find_last_of(CONTIG_SEPARATOR);
    if(colonIndex == -1) {
        contig = str;
        start = 1;
        end = INT32_MAX;
    } else {
        contig = str.substr(0, colonIndex);
        const int dashIndex = str.find(START_END_SEPARATOR, colonIndex);
        if(dashIndex == -1) {
            if(str.find_last_of(END_OF_CONTIG) == str.size() - END_OF_CONTIG.size()) {
                std::string pos = str.substr(colonIndex+1, str.size()-colonIndex-2);
                start = parsePosition(pos);
                end = INT32_MAX;
            } else {
                std::string pos = str.substr(colonIndex+1, str.size()-colonIndex-1);
                start = parsePosition(pos);
                end = start;
            }
        } else {
            std::string pos = str.substr(colonIndex+1, dashIndex-colonIndex-1);
            start = parsePosition(pos);
            pos = str.substr(dashIndex+1, str.size()-dashIndex-1);
            end = parsePosition(pos);
        }
    }

    validatePositions(contig, start, end);
    this->contig = contig;
    this->start = start;
    this->end = end;
}

void SimpleInterval::validatePositions(const std::string& contig, const int start, const int end) {
    if(contig.empty() || start <= 0 || start > end){
        throw std::invalid_argument("Argument input error.");
    }
}

void SimpleInterval::printfInterval() const {
    std::cout << "Interval" << start << "~" << end << std::endl;
}

bool SimpleInterval::isValid(const std::string& contig, const int start, const int end) {
    return (!contig.empty()) && start > 0 && end >= start;
}

int SimpleInterval::parsePosition(std::string pos) {
    int postion;
    pos = Mutect2Utils::replaceWith(pos, ",", "");
    std::stringstream ss;
    ss << pos;
    ss >> postion;
    if(ss.eof() && !ss.fail())
        return postion;
    else
        throw std::invalid_argument("Contig input error.");
}

bool SimpleInterval::operator==(const SimpleInterval &interval) const {
    if (contig == interval.contig && start == interval.start && end == interval.end)
        return true;
    else
        return false;
}

int SimpleInterval::hashCode() {
    std::hash<std::string> h;
    int result = start;
    result = 31 * result + end;
    result = 31 * result + h(contig);
    return result;
}

bool SimpleInterval::overlapsWithMargin(Locatable *other, const int margin) {
    Mutect2Utils::validateArg(margin >= 0, "Given margin is negative.");
    if( other == nullptr || other->getContig().empty())
        return false;
    else
        return (this->contig == other->getContig()) && this->start <= other->getEnd() + margin && other->getStart() - margin <= this->end;
}

bool SimpleInterval::overlaps(Locatable *other) {
    return overlapsWithMargin(other, 0);
}

SimpleInterval* SimpleInterval::intersect(Locatable *other) {
    Mutect2Utils::validateArg(overlaps(other), "SimpleInterval::intersect(): The two intervals need to overlap.");
    return new SimpleInterval(getContig(), std::max(start, other->getStart()), std::min(end, other->getEnd()));
}

SimpleInterval* SimpleInterval::mergeWithContiguous(Locatable* other){
    Mutect2Utils::validateArg(other != nullptr, "Null object is not allowed here.");
    Mutect2Utils::validateArg(contiguous(other), "The two intervals need to be contiguous.");
    return new SimpleInterval(getContig(), std::min(start, other->getStart()), std::max(end, other->getEnd()));
}

bool SimpleInterval::contiguous(Locatable *other) {
    Mutect2Utils::validateArg(other != nullptr, "Null object is not allowed here.");
    return contig == other->getContig() && start <= other->getEnd() + 1 && other->getStart() <= end + 1;
}

SimpleInterval* SimpleInterval::spanWith(Locatable *other) {
    Mutect2Utils::validateArg(other != nullptr, "Null object is not allowed here.");
    Mutect2Utils::validateArg(contig == other->getContig(), "Cannot get span for intervals on different contigs.");
    return new SimpleInterval(getContig(), std::min(start, other->getStart()), std::max(end, other->getEnd()));
}

SimpleInterval* SimpleInterval::expandWithinContig(const int padding, const int contigLength) {
    if(padding < 0)
        throw std::invalid_argument("Padding must be >= 0.");

    return IntervalUtils::trimIntervalToContig(contig, start - padding, end + padding, contigLength);
}