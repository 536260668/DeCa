//
// Created by 梦想家xixi on 2021/10/11.
//

#ifndef MUTECT2CPP_MASTER_SIMPLEINTERVAL_H
#define MUTECT2CPP_MASTER_SIMPLEINTERVAL_H

#endif //MUTECT2CPP_MASTER_SIMPLEINTERVAL_H

#include <string>

static const char CONTIG_SEPARATOR = ':';
static const char START_END_SEPARATOR = '-';
static const std::string END_OF_CONTIG = "+";

class SimpleInterval
{
private:
    static const long serialVersionUID = 1L;
    int start;
    int end;
    std::string contig;

public:
    SimpleInterval(int start, int end, std::string contig) : start(start), end(end), contig(contig) {}
    ~SimpleInterval() {}
    void printfInterval();
};