//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_SEQVERTEX_H
#define MUTECT2CPP_MASTER_SEQVERTEX_H

#include "BaseVertex.h"

class SeqVertex : public BaseVertex{
public:
    /**
     * Create a new SeqVertex with sequence and the next available id
     * @param sequence our base sequence
     */
    SeqVertex(uint8_t* sequence, int length) : BaseVertex(sequence, length) {}

    int getId() {return hashCode();}

    long hashCode() const;

    bool operator<(const SeqVertex & other) const;

    SeqVertex* withoutSuffix(uint8_t* suffix, int length);

    SeqVertex* withoutPrefixAndSuffix(uint8_t* prefix, int preLength, uint8_t* suffix, int sufLength);
};


#endif //MUTECT2CPP_MASTER_SEQVERTEX_H
