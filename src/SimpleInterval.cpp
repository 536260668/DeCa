//
// Created by 梦想家xixi on 2021/10/11.
//

#include "SimpleInterval.h"
#include <iostream>

SimpleInterval::SimpleInterval(const std::string &contig, int start, int end) : contig(ContigMap::getContigInt(contig)),
                                                                                start(start), end(end) {
	validatePositions(this->contig, start, end);
}

SimpleInterval::SimpleInterval(std::string &&contig, int start, int end) : contig(ContigMap::getContigInt(contig)),
                                                                           start(start), end(end) {
}

SimpleInterval::SimpleInterval(int contig, int start, int end) : contig(contig), start(start), end(end) {
}

SimpleInterval::SimpleInterval(SimpleInterval const &simpleInterval) : contig(simpleInterval.contig),
                                                                       start(simpleInterval.start),
                                                                       end(simpleInterval.end) {
	validatePositions(contig, start, end);
}

SimpleInterval::SimpleInterval(const std::shared_ptr<Locatable> &pLocatable) : contig(
		ContigMap::getContigInt(pLocatable->getContig())), start(pLocatable->getStart()), end(pLocatable->getEnd()) {}

SimpleInterval::SimpleInterval(std::string &str) {
	if (str.empty())
		throw std::invalid_argument("Null object is not allowed here.");

	std::string m_contig;
	int m_start;
	int m_end;

	const int colonIndex = str.find_last_of(CONTIG_SEPARATOR);
	if (colonIndex == -1) {
		m_contig = str;
		m_start = 1;
		m_end = INT32_MAX;
	} else {
		m_contig = str.substr(0, colonIndex);
		const int dashIndex = str.find(START_END_SEPARATOR, colonIndex);
		if (dashIndex == -1) {
			if (str.find_last_of(END_OF_CONTIG) == str.size() - END_OF_CONTIG.size()) {
				std::string pos = str.substr(colonIndex + 1, str.size() - colonIndex - 2);
				m_start = parsePosition(pos);
				m_end = INT32_MAX;
			} else {
				std::string pos = str.substr(colonIndex + 1, str.size() - colonIndex - 1);
				m_start = parsePosition(pos);
				m_end = m_start;
			}
		} else {
			std::string pos = str.substr(colonIndex + 1, dashIndex - colonIndex - 1);
			m_start = parsePosition(pos);
			pos = str.substr(dashIndex + 1, str.size() - dashIndex - 1);
			m_end = parsePosition(pos);
		}
	}

	int contigInt = ContigMap::getContigInt(m_contig);
	validatePositions(contigInt, m_start, m_end);
	contig = contigInt;
	start = m_start;
	end = m_end;
}

void SimpleInterval::clearContig() {
	this->contig = -1;
}

void SimpleInterval::validatePositions(int contig, const int start, const int end) {
	if (contig == -1 || start < 0 || start > end)
		throw std::invalid_argument("Argument input error.");
}

int SimpleInterval::parsePosition(std::string pos) {
	int postion;
	pos = Mutect2Utils::replaceWith(pos, ",", "");
	std::stringstream ss;
	ss << pos;
	ss >> postion;
	if (ss.eof() && !ss.fail())
		return postion;
	else
		throw std::invalid_argument("Contig input error.");
}

bool SimpleInterval::operator==(const SimpleInterval &interval) const {
	if (contig == interval.contig && start == interval.start && end == interval.end)
		return true;
	return false;
}

int SimpleInterval::hashCode() const {
	return 31 * 31 * start + 31 * end + contig;
}

bool SimpleInterval::overlapsWithMargin(const std::shared_ptr<Locatable> &other, const int margin) const {
	if (margin < 0)
		throw std::invalid_argument("Given margin is negative.");

	if (other == nullptr || other->getContig().empty())
		return false;

	return (this->contig == ContigMap::getContigInt(other->getContig())) && this->start <= other->getEnd() + margin &&
	       other->getStart() - margin <= this->end;
}

bool SimpleInterval::overlaps(const std::shared_ptr<Locatable> &other) {
	return overlapsWithMargin(other, 0);
}

std::shared_ptr<SimpleInterval> SimpleInterval::intersect(const std::shared_ptr<Locatable> &other) {
	if (!overlaps(other))
		throw std::invalid_argument("SimpleInterval::intersect(): The two intervals need to overlap.");

	return std::make_shared<SimpleInterval>(getContigInt(), std::max(start, other->getStart()),
	                                        std::min(end, other->getEnd()));
}

std::shared_ptr<SimpleInterval> SimpleInterval::mergeWithContiguous(const std::shared_ptr<Locatable> &other) {
	if (other == nullptr)
		throw std::invalid_argument("Null object is not allowed here.");
	if (!contiguous(other.get()))
		throw std::invalid_argument("The two intervals need to be contiguous.");

	return std::make_shared<SimpleInterval>(getContigInt(), std::min(start, other->getStart()),
	                                        std::max(end, other->getEnd()));
}

bool SimpleInterval::contiguous(Locatable *other) const {
	return contig == ContigMap::getContigInt(other->getContig()) && start <= other->getEnd() + 1 &&
	       other->getStart() <= end + 1;
}

std::shared_ptr<SimpleInterval> SimpleInterval::spanWith(const std::shared_ptr<Locatable> &other) const {
	if (other == nullptr)
		throw std::invalid_argument("Null object is not allowed here.");
	if (contig != ContigMap::getContigInt(other->getContig()))
		throw std::invalid_argument("Cannot get span for intervals on different contigs.");

	return std::make_shared<SimpleInterval>(getContigInt(), std::min(start, other->getStart()),
	                                        std::max(end, other->getEnd()));
}

std::shared_ptr<SimpleInterval> SimpleInterval::expandWithinContig(const int padding, const int contigLength) const {
	if (padding < 0)
		throw std::invalid_argument("Padding must be >= 0.");

	return IntervalUtils::trimIntervalToContig(ContigMap::getContigString(contig), start - padding, end + padding,
	                                           contigLength);
}

std::ostream &operator<<(std::ostream &os, const SimpleInterval &simpleInterval) {
	os << "contig:" << simpleInterval.contig << "  start:" << simpleInterval.start << "   end:" << simpleInterval.end
	   << std::endl;
	return os;
}


std::shared_ptr<SimpleInterval>
SimpleInterval::expandWithinContig(int padding, SAMSequenceDictionary *sequenceDictionary) const {
	if (sequenceDictionary == nullptr)
		throw std::invalid_argument("null is not allowed there");

	SAMSequenceRecord &contigRecord = sequenceDictionary->getSequences()[contig];
	return expandWithinContig(padding, contigRecord.getSequenceLength());
}

SimpleInterval::SimpleInterval() : contig(-1), start(0), end(0) {}

void SimpleInterval::printInfo() const {
	std::cout << getContig() << " " << getStart() + 1 << " " << getEnd() + 1 << std::endl;
}

