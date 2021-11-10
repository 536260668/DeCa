//
// Created by 梦想家xixi on 2021/10/12.
//

#ifndef MUTECT2CPP_MASTER_MUTECT2UTILS_H
#define MUTECT2CPP_MASTER_MUTECT2UTILS_H

#include <string>
#include <cfloat>
#include <vector>

static double POSITIVE_INFINITY = DBL_MAX;

class Mutect2Utils
{
public:
    static std::string replaceWith(std::string& str1, const std::string& str2, const std::string& str3);
    static bool overlaps(int start, int end, int start2, int end2);
    static bool encloses(int outerStart, int outerEnd, int innerStart, int innerEnd);
    static void validateArg(bool condition, std::string msg);
    static bool goodProbability(double result);
    static double logLikelihoodRatio(int nRef, std::vector<uint8_t> altQuals, int repeatFactor);
    static double logLikelihoodRatio(int refCount, int altCount, double errorProbability);
};

#endif //MUTECT2CPP_MASTER_MUTECT2UTILS_H