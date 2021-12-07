//
// Created by 梦想家xixi on 2021/11/8.
//

#include <cstring>
#include "Allele.h"
#include <stdexcept>
#include "StringUtils.h"
const uint8_t* Allele::EMPTY_ALLELE_BASES = new uint8_t[0];
const std::string Allele::NO_CALL_STRING = ".";
const std::string Allele::SPAN_DEL_STRING = "*";
const std::string Allele::NON_REF_STRING = "<NON_REF>";
const std::string Allele::UNSPECIFIED_ALTERNATE_ALLELE_STRING = "<*>";
Allele Allele::REF_A(new uint8_t[1]{'A'}, 1, true);
Allele Allele::ALT_A(new uint8_t[1]{'A'}, 1, false);
Allele Allele::REF_C(new uint8_t[1]{'C'}, 1, true);
Allele Allele::ALT_C(new uint8_t[1]{'C'}, 1, false);
Allele Allele::REF_G(new uint8_t[1]{'G'}, 1, true);
Allele Allele::ALT_G(new uint8_t[1]{'G'}, 1, false);
Allele Allele::REF_T(new uint8_t[1]{'T'}, 1, true);
Allele Allele::ALT_T(new uint8_t[1]{'T'}, 1, false);
Allele Allele::REF_N(new uint8_t[1]{'N'}, 1, true);
Allele Allele::ALT_N(new uint8_t[1]{'N'}, 1, false);
Allele Allele::SPAN_DEL(new uint8_t[1]{'*'}, 1, false);
Allele Allele::NO_CALL(new uint8_t[1]{'.'}, 1, false);
Allele Allele::NON_REF_ALLELE(new uint8_t[9]{'<', 'N', 'O', 'N', '_', 'R', 'E', 'F', '>'}, 9, false);
Allele Allele::UNSPECIFIED_ALTERNATE_ALLELE(new uint8_t[3]{'<', '*', '>'}, 3, true);

Allele::Allele(uint8_t *bases, int length, bool isRef) : isRef(false), bases(nullptr), length(0){
    if(wouldBeNullAllele(bases, length)) {
        throw std::invalid_argument("Null alleles are not supported");
    } else if(wouldBeNoCallAllele(bases, length)) {
        this->bases = const_cast<uint8_t *>(EMPTY_ALLELE_BASES);
        this->isNoCall = true;
        if(isRef) {
            throw std::invalid_argument("Cannot tag a NoCall allele as the reference allele");
        } else {
            if(wouldBeSymbolicAllele(bases, length)) {
                this->isSymbolic = true;
                if(isRef) {
                    throw std::invalid_argument("Cannot tag a symbolic allele as the reference allele");
                } else {
                    StringUtils::toUpperCase(bases, length);
                }

                this ->isRef = isRef;
                this->bases = bases;
                if(!acceptableAlleleBases(bases, length, isRef)) {
                    throw std::invalid_argument("Unexpected base in allele bases");
                }
            }
        }
    }
}

bool Allele::wouldBeNullAllele(const uint8_t *bases, int length){
    return length == 1 && bases[0] == 45 || length == 0;
}

bool Allele::wouldBeNoCallAllele(const uint8_t *bases, int length) {
    return length == 1 && bases[0] == 46;
}

bool Allele::wouldBeSymbolicAllele(const uint8_t *bases, int length) {
    if(length <= 1) {
        return false;
    } else {
        return bases[0] == 60 || bases[length - 1] == 62 || wouldBeBreakpoint(bases, length) || wouldBeSingleBreakend(bases, length);
    }
}

bool Allele::wouldBeBreakpoint(const uint8_t *bases, int length) {
    if (length <= 1) {
        return false;
    } else {
        for(int i = 0; i < length; ++i) {
            uint8_t base = bases[i];
            if (base == 93 || base == 91) {
                return true;
            }
        }
        return false;
    }
}

bool Allele::wouldBeSingleBreakend(const uint8_t *bases, int length) {
    if (length <= 1) {
        return false;
    } else {
        return bases[0] == 46 || bases[length - 1] == 46;
    }
}

bool Allele::acceptableAlleleBases(const uint8_t *bases, int length, bool isReferenceAllele) {
    if (wouldBeNullAllele(bases, length)) {
        return false;
    } else if (!wouldBeNoCallAllele(bases, length) && !wouldBeSymbolicAllele(bases, length)) {
        if (wouldBeStarAllele(bases, length)) {
            return !isReferenceAllele;
        } else {
            uint8_t* var2 = const_cast<uint8_t *>(bases);
            int var3 = length;
            int var4 = 0;

            while(var4 < var3) {
               uint8_t base = var2[var4];
                switch(base) {
                    case 65:
                    case 67:
                    case 71:
                    case 78:
                    case 84:
                    case 97:
                    case 99:
                    case 103:
                    case 110:
                    case 116:
                        ++var4;
                        break;
                    default:
                        return false;
                }
            }

            return true;
        }
    } else {
        return true;
    }
}

bool Allele::wouldBeStarAllele(const uint8_t *bases, int length) {
    return length == 1 && bases[0] == 42;
}

Allele::Allele(Allele &allele, bool ignoreRefState) : bases(allele.bases), isNoCall(allele.isNoCall), isSymbolic(allele.isSymbolic){
    this->isRef = !ignoreRefState && allele.isRef;
}

Allele *Allele::create(uint8_t *bases, int length, bool isRef) {
    if(bases == nullptr) {
        throw std::invalid_argument("create: the Allele base string cannot be null; use new Allele() or new Allele(\"\") to create a Null allele");
    } else if (length == 1) {
        switch (bases[0]) {
            case 42:
                if(isRef) {
                    throw std::invalid_argument("Cannot tag a spanning deletions allele as the reference allele");
                }
                return &SPAN_DEL;

            case 46:
                if(isRef) {
                    throw std::invalid_argument("Cannot tag a NoCall allele as the reference allele");
                }
                return &NO_CALL;

            case 65:
            case 97:
                return isRef ? &REF_A : &ALT_A;
            case 67:
            case 99:
                return isRef ? &REF_C : &ALT_C;
            case 71:
            case 103:
                return isRef ? &REF_G : &ALT_G;
            case 78:
            case 110:
                return isRef ? &REF_N : &ALT_N;
            case 84:
            case 116:
                return isRef ? &REF_T : &ALT_T;
            default:
                throw std::invalid_argument("Illegal base seen in the allele");
        }
    } else {
        return new Allele(bases, length, isRef);
    }
}

Allele *Allele::create(uint8_t base, bool isRef) {
    uint8_t * bases = new uint8_t[1]{base};
    Allele* res = create(bases, 1, isRef);
    delete[] bases;
    return res;
}

Allele *Allele::create(uint8_t base) {
    return create(base, false);
}

Allele *Allele::extend(Allele *left, uint8_t *right, int length) {
    if(left->isSymbolic) {
        throw std::invalid_argument("Cannot extend a symbolic allele");
    } else {
        uint8_t * bases = new uint8_t[left->length + length];
        memcpy(bases, left->bases, left->length);
        memcpy(bases+left->length, right, length);
        return create(bases, left->length+length, left->isRef);
    }
}

bool Allele::operator<(const Allele &other) const {
    if(isRef < other.getIsReference())
        return true;
    if(isRef > other.getIsReference())
        return false;
    if(isNoCall < other.getIsNoCall())
        return true;
    if(isNoCall < other.getIsNoCall())
        return false;
    for(int i = 0; i < length; i++) {
        if(bases[i] < other.bases[i])
            return true;
        if(bases[i] > other.bases[i])
            return false;
    }
    return true;

}

int Allele::getLength() const {
    return isSymbolic ? 0 : length;
}

bool Allele::operator==(const Allele &other) const {
    if(this == &other)
        return true;
    if(this->isRef != other.getIsReference() || this->isNoCall != other.getIsNoCall())
        return false;
    if(this->bases == other.getBases())
        return true;
    if(this->getLength() != other.getLength())
        return false;
    for(int i = 0; i < length; i++) {
        if(bases[i] != other.bases[i])
            return false;
    }
    return true;
}

std::string Allele::getBaseString() {
    char* tmp = new char[length+1]{0};
    memcpy(tmp, bases, length);
    std::string ret = tmp;
    delete[] tmp;
    return ret;
}

bool Allele::equals(Allele &other, bool ignoreRefState) {
    if(this == &other)
        return true;
    if(isRef != other.isRef && !ignoreRefState)
        return false;
    if(isNoCall != other.isNoCall)
        return false;
    if(bases == other.bases)
        return true;
    else{
        if(length != other.length)
            return false;
        for(int i = 0; i < length; i++) {
            if(bases[i] != other.bases[i]) {
                return false;
            }
        }
        return true;
    }
}
