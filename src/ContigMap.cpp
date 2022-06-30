//
// Created by hlf on 6/29/22.
//
#include "ContigMap.h"
#include "boost/utility.hpp"

std::vector<std::string> ContigMap::intToString;
std::unordered_map<std::string, int, hash_contig, equal_contig> ContigMap::stringToInt;
int ContigMap::mapSize;

void ContigMap::initial(int reserveSize) {
	mapSize = reserveSize;
	intToString.resize(mapSize + 1);
	stringToInt.reserve(mapSize + 1);
	stringToInt.insert(std::make_pair("", -1));
}

int ContigMap::getContigInt(const std::string &key) {
	auto iter = stringToInt.find(key);
	if (BOOST_UNLIKELY(iter == stringToInt.end()))
		throw std::invalid_argument("ContigMap key " + key + " not found.");
	return iter->second;
}

std::string ContigMap::getContigString(int key) {
	if (BOOST_UNLIKELY(key >= mapSize || key < -1))
		throw std::invalid_argument("ContigMap key " + std::to_string(key) + " not found.");
	return intToString[key + 1];
}

void ContigMap::insertPair(int intKey, const std::string &stringKey) {
	if (BOOST_UNLIKELY(intKey >= mapSize || intKey < -1))
		throw std::invalid_argument("ContigMap key " + std::to_string(intKey) + " not found.");
	if (BOOST_UNLIKELY(stringToInt.find(stringKey) != stringToInt.end()))
		throw std::invalid_argument("can't insert the same int-string pair twice.");
	intToString[intKey + 1] = stringKey;
	stringToInt.insert(std::make_pair(stringKey, intKey));
}

bool equal_contig::operator()(const std::string &s1, const std::string &s2) const {
	if (s1.length() != s2.length())
		return false;
	for (int i = 0; i < s1.length(); ++i) {
		if (s1[i] != s2[i])
			return false;
	}
	return true;
}

xxh::hash64_t hash_contig::operator()(const std::string &s1) const {
	return xxh::xxhash3<64>(s1);
}
