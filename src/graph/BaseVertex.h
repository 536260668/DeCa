//
// Created by 梦想家xixi on 2021/10/20.
//

#ifndef MUTECT2CPP_MASTER_BASEVERTEX_H
#define MUTECT2CPP_MASTER_BASEVERTEX_H


#include <cstdint>
#include <string>
#include <iostream>

class BaseVertex {
private:
    int cashedHashCode;
    int length;
    std::string additionalInfo;

    static int hashCode(uint8_t * a, int length);

protected:
    uint8_t* sequence;

public:
    /**
     * Create a new sequence vertex with sequence
     *
     * This code doesn't copy sequence for efficiency reasons, so sequence must absolutely not be modified
     * in any way after passing this sequence to the BaseVertex
     *
     * @param sequence a non-null sequence of bases contained in this vertex
     */
    BaseVertex(uint8_t* sequence, int length);

    virtual ~BaseVertex() = default;

    /**
     * Does this vertex have an empty sequence?
     *
     * That is, is it a dummy node that's only present for structural reasons but doesn't actually
     * contribute to the sequence of the graph?
     *
     * @return true if sequence is empty, false otherwise
     */
    bool isEmpty() const;

    /**
     * Get the length of this sequence
     * @return a positive integer >= 1
     */
    int getLength() const {return length;}

    bool operator==(const BaseVertex & other) const;

    bool operator<(const BaseVertex & other) const;

    int getHashCode() const {return cashedHashCode;}

    friend  std::ostream & operator<<(std::ostream &os, const BaseVertex & baseVertex);

    uint8_t * getSequence() const {return sequence;}

    void setAdditionalInfo(const std::string &info);

    std::string getAdditionalInfo() const{return additionalInfo;}

    bool hasAmbiguousSequence();

    bool seqEquals(BaseVertex* other);

    virtual uint8_t * getAdditionalSequence(bool source) {return getSequence();}
};


#endif //MUTECT2CPP_MASTER_BASEVERTEX_H
