//
// Created by hlf on 6/29/22.
//

#ifndef MUTECT2CPP_MASTER_CONTIGMAP_H
#define MUTECT2CPP_MASTER_CONTIGMAP_H

#include <unordered_map>
#include <string>
#include <stdexcept>

class ContigMap {
private:
	static std::unordered_map<int, std::string> intToString;
	static std::unordered_map<std::string, int> stringToInt;

public:
	static void initial(int reserveSize);

	static int getContigInt(const std::string &key);

	static std::string getContigString(int key);

	static void insertPair(int intKey, const std::string &stringKey);
};

#endif //MUTECT2CPP_MASTER_CONTIGMAP_H
