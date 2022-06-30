//
// Created by hlf on 6/29/22.
//
#include "ContigMap.h"

std::unordered_map<int, std::string> ContigMap::intToString;
std::unordered_map<std::string, int> ContigMap::stringToInt;

void ContigMap::initial(int reserveSize) {
	intToString.reserve(reserveSize + 1);
	stringToInt.reserve(reserveSize + 1);
	insertPair(-1, "");
}

int ContigMap::getContigInt(const std::string &key) {
	if (key.empty())
		return -1;
	if (stringToInt.find(key) == stringToInt.end())
		throw std::invalid_argument("ContigMap key " + key + " not found.");
	return stringToInt[key];
}

std::string ContigMap::getContigString(int key) {
	if (key == -1)
		return "";
	if (intToString.find(key) == intToString.end())
		throw std::invalid_argument("ContigMap key " + std::to_string(key) + " not found.");
	return intToString[key];
}

void ContigMap::insertPair(int intKey, const std::string &stringKey) {
	if (stringToInt.find(stringKey) != stringToInt.end() || intToString.find(intKey) != intToString.end())
		throw std::invalid_argument("can't insert the same int-string pair twice.");
	intToString.insert(std::make_pair(intKey, stringKey));
	stringToInt.insert(std::make_pair(stringKey, intKey));
}
