//
// Created by 梦想家xixi on 2021/11/12.
//

#ifndef MUTECT2CPP_MASTER_SMITHWATERMANALIGNER_H
#define MUTECT2CPP_MASTER_SMITHWATERMANALIGNER_H

#include <cstdint>
#include "SWParameters.h"
#include "SmithWatermanAlignment.h"

enum SWOverhangStrategy{SOFTCLIP, INDEL, LEADING_INDEL, IGNORE};

class SmithWatermanAligner {
public:
    static SWParameters ORIGINAL_DEFAULT;
    static SWParameters STANDARD_NGS;

    virtual SmithWatermanAlignment* align(std::shared_ptr<uint8_t[]> ref, int refLength, std::shared_ptr<uint8_t[]> alt, int altLength, SWParameters* parameters, SWOverhangStrategy overhangStrategy) = 0;
};


#endif //MUTECT2CPP_MASTER_SMITHWATERMANALIGNER_H
