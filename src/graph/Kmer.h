//
// Created by 梦想家xixi on 2021/10/19.
//

#ifndef MUTECT2CPP_MASTER_KMER_H
#define MUTECT2CPP_MASTER_KMER_H


#include <cstdint>

class Kmer {
private:
    uint8_t  *bases;
    int start;
    int length;
    int hash;

    /**
     *  Compute the hashcode for a KMer.
     *  Equivalent to <code>new String(bases, start, length).hashCode()</code>
     */
    static int hashCode(const uint8_t* bases, int start, int length);

public:
    Kmer(uint8_t* kmer, int length);

    Kmer(Kmer const & kmer);

    /**
     * Create a new kmer backed by the bases in bases, spanning start -> start + length
     *
     * Under no circumstances can bases be modified anywhere in the client code.  This does not make a copy
     * of bases for performance reasons
     *
     * @param bases an array of bases
     * @param start the start of the kmer in bases, must be >= 0 and < bases.length
     * @param length the length of the kmer.  Must be >= 0 and start + length < bases.length
     */
    Kmer(uint8_t* kmer, int start, int length);

    Kmer subKmer(int newStart, int newLength);

    //delete
    uint8_t * getBases() const;

    /**
     * The length of this kmer
     * @return an integer >= 0
     */
    int getLength() const {return length;}

    int getDifferingPositions(Kmer other, int maxDistance, int * differingIndeces, uint8_t * differingBases);

    bool operator<(const Kmer & other) const;

    bool operator==(const Kmer & other) const;
};


#endif //MUTECT2CPP_MASTER_KMER_H
