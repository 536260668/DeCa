//
// Created by lhh on 11/13/21.
//

#ifndef MUTECT2CPP_MASTER_COMMONINFO_H
#define MUTECT2CPP_MASTER_COMMONINFO_H

#include <string>
#include <set>

/**
 * Common utility routines for VariantContext and Genotype
 */
class CommonInfo {
public:
    constexpr static double NO_LOG10_PERROR = 1.0;

    CommonInfo(std::string name, double log10PError, std::set<std::string> filters);

    void setLog10PError(double log10PError);


private:
    double log10PError = NO_LOG10_PERROR;
    std::string name;
    std::set<std::string> filters;
    // map<string, object> attributes;  // TODO: add attributes field
};


#endif //MUTECT2CPP_MASTER_COMMONINFO_H
