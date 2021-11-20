//
// Created by 梦想家xixi on 2021/11/8.
//

#ifndef MUTECT2CPP_MASTER_ALLELE_H
#define MUTECT2CPP_MASTER_ALLELE_H


#include <cstdint>
#include <string>

class Allele {
private:
    static const uint8_t* EMPTY_ALLELE_BASES;
    static const char SINGLE_BREAKEND_INDICATOR = '.';
    static const char BREAKEND_EXTENDING_RIGHT = '[';
    static const char BREAKEND_EXTENDING_LEFT = ']';
    static const char SYMBOLIC_ALLELE_START = '<';
    static const char SYMBOLIC_ALLELE_END = '>';
    bool isRef;
    bool isNoCall;
    bool isSymbolic;
    uint8_t *bases;
    int length;

    static Allele REF_A;
    static Allele ALT_A;
    static Allele REF_C;
    static Allele ALT_C;
    static Allele REF_G;
    static Allele ALT_G;
    static Allele REF_T;
    static Allele ALT_T;
    static Allele REF_N;
    static Allele ALT_N;
    static Allele SPAN_DEL;
    static Allele NO_CALL;
    static Allele NON_REF_ALLELE;
    static Allele UNSPECIFIED_ALTERNATE_ALLELE;
    static Allele SV_SIMPLE_DEL;
    static Allele SV_SIMPLE_INS;
    static Allele SV_SIMPLE_INV;
    static Allele SV_SIMPLE_CNV;
    static Allele SV_SIMPLE_DUP;


public:
    static const std::string NO_CALL_STRING;
    static const std::string SPAN_DEL_STRING;
    static const std::string NON_REF_STRING;
    static const std::string UNSPECIFIED_ALTERNATE_ALLELE_STRING;
    static bool wouldBeNullAllele(const uint8_t * bases, int length);
    static bool wouldBeNoCallAllele(const uint8_t * bases, int length);
    static bool wouldBeSymbolicAllele(const uint8_t* bases, int length);
    static bool wouldBeBreakpoint(const uint8_t* bases, int length);
    static bool wouldBeSingleBreakend(const uint8_t* bases, int length);
    static bool wouldBeStarAllele(const uint8_t* bases, int length);
    static bool acceptableAlleleBases(const uint8_t* bases, int length, bool isReferenceAllele);
    static Allele* create(uint8_t * bases, int length, bool isRef);
    static Allele* create(uint8_t base, bool isRef);
    static Allele* create(uint8_t base);
    static Allele* extend(Allele * left, uint8_t * right, int length);
    bool getIsNoCall() const  {return isNoCall;}
    bool getIsCalled() const {return !isNoCall;}
    bool getIsReference() const {return isRef;}
    bool getIsNonReference() const {return !isRef;}
    bool getIsSymbolic() const {return isSymbolic;}
    bool getIsBreakpoint() const {return wouldBeBreakpoint(bases, length);}
    bool getIsSingleBreakend() const {return wouldBeSingleBreakend(bases, length);}
    uint8_t* getBases() const {return bases;}
    bool operator<(const Allele & other) const;
    int getLength() const;
    int getBasesLength() const {return length;}

protected:
    Allele(uint8_t* bases, int length, bool isRef);
    Allele(Allele &  allele, bool ignoreRefState);
};


#endif //MUTECT2CPP_MASTER_ALLELE_H
