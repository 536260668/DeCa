//
// Created by 梦想家xixi on 2021/10/11.
//

#include "SimpleInterval.h"
#include <iostream>

SimpleInterval::SimpleInterval(const std::string& contig, int start, int end) : contig(contig), start(start), end(end){
    validatePositions(contig, start, end);
}

SimpleInterval::SimpleInterval(std::string&& contig, int start, int end) : contig(contig), start(start), end(end){

}

SimpleInterval::SimpleInterval(SimpleInterval const &simpleInterval) : contig(simpleInterval.contig), start(simpleInterval.start), end(simpleInterval.end){
    validatePositions(contig, start, end);
}

SimpleInterval::SimpleInterval(std::string& str) {
    Mutect2Utils::validateArg(!str.empty(), "Null object is not allowed here.");

    std::string m_contig;
    int m_start;
    int m_end;

    const int colonIndex = str.find_last_of(CONTIG_SEPARATOR);
    if(colonIndex == -1) {
        m_contig = str;
        m_start = 1;
        m_end = INT32_MAX;
    } else {
        m_contig = str.substr(0, colonIndex);
        const int dashIndex = str.find(START_END_SEPARATOR, colonIndex);
        if(dashIndex == -1) {
            if(str.find_last_of(END_OF_CONTIG) == str.size() - END_OF_CONTIG.size()) {
                std::string pos = str.substr(colonIndex+1, str.size()-colonIndex-2);
                m_start = parsePosition(pos);
                m_end = INT32_MAX;
            } else {
                std::string pos = str.substr(colonIndex+1, str.size()-colonIndex-1);
                m_start = parsePosition(pos);
                m_end = m_start;
            }
        } else {
            std::string pos = str.substr(colonIndex+1, dashIndex-colonIndex-1);
            m_start = parsePosition(pos);
            pos = str.substr(dashIndex+1, str.size()-dashIndex-1);
            m_end = parsePosition(pos);
        }
    }

    validatePositions(m_contig, m_start, m_end);
    contig = m_contig;
    start = m_start;
    end = m_end;
}

void SimpleInterval::clearContig()
{
    this->contig.clear();
}

void SimpleInterval::validatePositions(const std::string& contig, const int start, const int end) {
    if(contig.empty() || start < 0 || start > end){
        throw std::invalid_argument("Argument input error.");
    }
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

int SimpleInterval::hashCode() const{
    std::hash<std::string> h;
    int result = start;
    result = 31 * result + end;
    result = 31 * result + h(contig);
    return result;
}

bool SimpleInterval::overlapsWithMargin(const std::shared_ptr<Locatable> & other, const int margin) const {
    Mutect2Utils::validateArg(margin >= 0, "Given margin is negative.");
    if( other == nullptr || other->getContig().empty())
        return false;
    else
        return (this->contig == other->getContig()) && this->start <= other->getEnd() + margin && other->getStart() - margin <= this->end;
}

bool SimpleInterval::overlaps(const std::shared_ptr<Locatable> & other) {
    return overlapsWithMargin(other, 0);
}

std::shared_ptr<SimpleInterval> SimpleInterval::intersect(const std::shared_ptr<Locatable> & other) {
    Mutect2Utils::validateArg(overlaps(other), "SimpleInterval::intersect(): The two intervals need to overlap.");
    std::shared_ptr<SimpleInterval> ret = std::make_shared<SimpleInterval>(getContig(), std::max(start, other->getStart()), std::min(end, other->getEnd()));
    return ret;
}

std::shared_ptr<SimpleInterval> SimpleInterval::mergeWithContiguous(const std::shared_ptr<Locatable> & other){
    Mutect2Utils::validateArg(other != nullptr, "Null object is not allowed here.");
    Mutect2Utils::validateArg(contiguous(other.get()), "The two intervals need to be contiguous.");
    std::shared_ptr<SimpleInterval> ret = std::make_shared<SimpleInterval>(getContig(), std::min(start, other->getStart()), std::max(end, other->getEnd()));
    return ret;
}

bool SimpleInterval::contiguous(Locatable *other) {
    Mutect2Utils::validateArg(other != nullptr, "Null object is not allowed here.");
    return contig == other->getContig() && start <= other->getEnd() + 1 && other->getStart() <= end + 1;
}

std::shared_ptr<SimpleInterval> SimpleInterval::spanWith(const std::shared_ptr<Locatable> &other) {
    Mutect2Utils::validateArg(other != nullptr, "Null object is not allowed here.");
    Mutect2Utils::validateArg(contig == other->getContig(), "Cannot get span for intervals on different contigs.");
    std::shared_ptr<SimpleInterval> ret = std::make_shared<SimpleInterval>(getContig(), std::min(start, other->getStart()), std::max(end, other->getEnd()));
    return ret;
}

std::shared_ptr<SimpleInterval> SimpleInterval::expandWithinContig(const int padding, const int contigLength) {
    if(padding < 0)
        throw std::invalid_argument("Padding must be >= 0.");

    return IntervalUtils::trimIntervalToContig(contig, start - padding, end + padding, contigLength);
}

std::ostream & operator<<(std::ostream &os, const SimpleInterval& simpleInterval) {
    os << "contig:" << simpleInterval.contig << "  start:" << simpleInterval.start << "   end:" << simpleInterval.end << std::endl;
    return os;
}

SimpleInterval::SimpleInterval(const std::shared_ptr<Locatable> & pLocatable) : contig(pLocatable->getContig()), start(pLocatable->getStart()), end(pLocatable->getEnd()){}

std::shared_ptr<SimpleInterval> SimpleInterval::expandWithinContig(int padding, SAMSequenceDictionary *sequenceDictionary) {
    Mutect2Utils::validateArg(sequenceDictionary, "null is not allowed there");
    SAMSequenceRecord& contigRecord = sequenceDictionary->getSequence(contig);
    return expandWithinContig(padding, contigRecord.getSequenceLength());
}

SimpleInterval::SimpleInterval() : start(0), end(0), contig(""){

}

