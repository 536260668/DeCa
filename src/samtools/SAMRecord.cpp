//
// Created by 梦想家xixi on 2021/12/20.
//

#include "SAMRecord.h"
#include "SAMFlag.h"
#include "ReadConstants.h"
#include "SAMUtils.h"
#include "ReadUtils.h"
#include <sstream>
#include <iostream>

const std::string SAMRecord::NO_ALIGNMENT_REFERENCE_NAME = "*";

std::string &SAMRecord::getName() {
    return mReadName;
}

void SAMRecord::setName(std::string& name) {
    mReadName = name;
}

int SAMRecord::getLength() const {
    return baseLength;
}

void SAMRecord::setPosition(std::string& contig, int start) {
    if(contig.empty() || contig == NO_ALIGNMENT_REFERENCE_NAME || start < 1) {
        throw std::invalid_argument("contig must be non-null and start must be >= 1");
    }
    mReferenceName = contig;
    mAlignmentStart = start;
    mAlignmentEnd = start + mCigar->getReferenceLength();
    setFlag(false, 4);
}

void SAMRecord::setFlag(bool flag, int bit) {
    if(flag) {
        mFlags |= bit;
    } else {
        mFlags &= ~bit;
    }
}

void SAMRecord::setPosition(Locatable *locatable) {
    Mutect2Utils::validateArg(locatable, "Cannot set read position to null");
    std::string contig = locatable->getContig();
    setPosition(contig, locatable->getStart());
}

std::string &SAMRecord::getAssignedContig() {
    return mReferenceName;
}

int SAMRecord::getAssignedStart() {
    return mAlignmentStart;
}

bool SAMRecord::getReadUnmappedFlag() {
    return (mFlags & 4) != 0;
}

bool SAMRecord::isUnmapped() {
    return getReadUnmappedFlag() || mReferenceName.empty() || mReferenceName == NO_ALIGNMENT_REFERENCE_NAME ||
    mAlignmentStart == NO_ALIGNMENT_START;
}

int SAMRecord::getUnclippedStart() {
    if(isUnmapped()) {
        return ReadConstants::UNSET_POSITION;
    }
    return SAMUtils::getUnclippedStart(mAlignmentStart, mCigar);
}

int SAMRecord::getUnclippedEnd() {
    if(isUnmapped())
        return ReadConstants::UNSET_POSITION;

    return SAMUtils::getUnclippedEnd(mAlignmentEnd, mCigar);
}

std::string& SAMRecord::getMateContig() {
    if(isUnmapped())
        return (std::string&)"";
    return mMateReferenceName;
}

bool SAMRecord::isPaired() {
    return (mFlags & 1) != 0;
}

int SAMRecord::getMateStart() {
    if(mateIsUnmapped()) {
        return ReadConstants::UNSET_POSITION;
    }
    return mMateAlignmentStart;
}

bool SAMRecord::mateIsUnmapped() {
    Mutect2Utils::validateArg(isPaired(), "Cannot get mate information for an unpaired read");

    return getMateUnmappedFlag() || mMateReferenceName.empty() || mMateReferenceName == NO_ALIGNMENT_REFERENCE_NAME
    || mMateAlignmentStart == NO_ALIGNMENT_START;
}

void SAMRecord::requireReadPaired() {
    if((mFlags & 1) == 0) {
        throw std::invalid_argument("Inappropriate call if not paired read");
    }
}

bool SAMRecord::getMateUnmappedFlagUnchecked() {
    return (mFlags & 8) != 0;
}

bool SAMRecord::getMateUnmappedFlag() {
    requireReadPaired();
    return getMateUnmappedFlagUnchecked();
}

void SAMRecord::setMatePosition(std::string &contig, int start) {
    if(!contig.empty() || contig == NO_ALIGNMENT_REFERENCE_NAME || start < 1) {
        throw std::invalid_argument("contig must be non-null and start must be >= 1");
    }

    setIsPaired(true);
    //TODO::实现setMateReferenceName
    mMateReferenceName = contig;
    mMateAlignmentStart = start;
    setFlag(false, 8);
}

void SAMRecord::setReadPairedFlag(bool flag) {
    setFlag(flag, 1);
}

void SAMRecord::setProperPairFlag(bool flag) {
    setFlag(flag, 2);
}

void SAMRecord::setIsPaired(bool isPaired) {
    setReadPairedFlag(isPaired);
    if(! isPaired) {
        setProperPairFlag(false);
    }
}

void SAMRecord::setMatePosition(Locatable *locatable) {
    Mutect2Utils::validateArg(locatable, "Cannot set mate position to null");
    std::string contig = locatable->getContig();
    setMatePosition(contig, locatable->getStart());
}

int SAMRecord::getFragmentLength() const {
    return mInferredInsertSize;
}

void SAMRecord::setFragmentLength(int fragmentLength) {
    mInferredInsertSize = fragmentLength;
}

int SAMRecord::getMappingQuality() const {
    return mMappingQuality != NO_MAPPING_QUALITY ? mMappingQuality : NO_MAPPING_QUALITY;
}

void SAMRecord::setMappingQuality(int mappingQuality) {
    Mutect2Utils::validateArg(mappingQuality >= 0 && mappingQuality <= 255, "mapping quality must be >= 0 and <= 255");
    mMappingQuality = mappingQuality;
}

uint8_t *SAMRecord::getBases() {
    if(mReadBases != nullptr){
        uint8_t * ret = new uint8_t[baseLength+1]{0};
        std::copy(mReadBases, mReadBases + baseLength, ret);
        return ret;
    } else {
        return nullptr;
    }
}

uint8_t *SAMRecord::getBasesNoCopy() {
    return mReadBases;
}

void SAMRecord::setBases(uint8_t *bases, int length) {
    delete[] bases;
    mReadBases = bases;
    baseLength = length;
}

uint8_t *SAMRecord::getBaseQualities() {
    if(mBaseQualities != nullptr){
        uint8_t * ret = new uint8_t[baseQualitiesLength+1]{0};
        std::copy(mBaseQualities, mBaseQualities+baseQualitiesLength, ret);
        return ret;
    } else {
        return nullptr;
    }
}

uint8_t *SAMRecord::getBaseQualitiesNoCopy() {
    return mBaseQualities;
}

int SAMRecord::getLength() {
    return baseLength;
}

int SAMRecord::getBaseQualitiesLength() {
    return baseQualitiesLength;
}

void SAMRecord::setBaseQualities(uint8_t *baseQualities, int length) {
    delete[] mBaseQualities;
    mBaseQualities = baseQualities;
    baseQualitiesLength = length;
}

Cigar *SAMRecord::getCigar() {
    //TODO:验证是否需要返回拷贝后的cigar
    return mCigar;
}

std::vector<CigarElement>& SAMRecord::getCigarElements() {
    return mCigar->getCigarElements();
}

CigarElement SAMRecord::getCigarElement(const int index) {
    return mCigar->getCigarElement(index);
}

void SAMRecord::setCigar(Cigar *cigar) {
    //TODO::实现AlignmentBlock
    isCalAdaptorBoundary = false;
    delete mCigar;
    mCigar = cigar;
    mAlignmentEnd = mAlignmentStart + mCigar->getReferenceLength();
}

bool SAMRecord::getProperPairFlag() {
    requireReadPaired();
    return getProperPairFlagUnchecked();
}

bool SAMRecord::getProperPairFlagUnchecked() const {
    return (mFlags & 2) != 0;
}

bool SAMRecord::isProperlyPaired() {
    return isPaired() && getProperPairFlag();
}

void SAMRecord::setIsProperlyPaired(bool isProperlyPaired) {
    if(isProperlyPaired) {
        setIsPaired(true);
    }

    setProperPairFlag(isProperlyPaired);
}

SAMRecord::SAMRecord(uint8_t *base, int baseLength, uint8_t *baseQualities, int baseQualitiesLength,
                     std::string &name) : mReadBases(base), baseLength(baseLength), mBaseQualities(baseQualities),baseQualitiesLength(baseQualitiesLength), mReadName(name){}

int SAMRecord::getStart() {
    if(isUnmapped()) {
        return ReadConstants::UNSET_POSITION;
    }
    return mAlignmentStart;
}

int SAMRecord::getEnd() {
    if(isUnmapped()) {
        return ReadConstants::UNSET_POSITION;
    }

    return mAlignmentEnd;
}

void SAMRecord::setIsUnmapped() {
    setReadUnmappedFlag(true);
}

void SAMRecord::setReadUnmappedFlag(bool flag) {
    setFlag(flag, 4);
}

void SAMRecord::clearAttributes() {
    mAttributes = nullptr;
}

void SAMRecord::setAttribute(short tag, void *value, Void_Type type, int length, bool isUnsignedArray) {
    if(value == nullptr) {
        if(mAttributes != nullptr) {
            mAttributes = SAMBinaryTagAndValue::remove(mAttributes, tag);
        }
    } else {
        SAMBinaryTagAndValue* tmp;
        if(!isUnsignedArray) {
            tmp = new SAMBinaryTagAndValue(tag, value, type, length);
        } else {
            if(type == Short_Array_Type || type == Int_Array_Type || type == Float_Array_Type || type == Uint8_t_Array_Type) {
                throw std::invalid_argument("Attribute type cannot be encoded as an unsigned array.");
            } else {
                tmp = new SAMBinaryTagAndValue(tag, value, type, length);
            }
            if(mAttributes == nullptr) {
                mAttributes == tmp;
            } else {
                mAttributes = SAMBinaryTagAndValue::insert(mAttributes, tmp);
            }
        }

    }
}

void *SAMRecord::getAttribute(short tag) {
    if(mAttributes == nullptr) {
        return nullptr;
    } else {
        SAMBinaryTagAndValue* tmp = mAttributes->find(tag);
        return tmp != nullptr ? tmp->value : nullptr;
    }
}

std::string &SAMRecord::getReadGroup() {
    return *(std::string*)getAttribute(SAMUtils::makeBinaryTag((std::string &) "RG"));
}

void SAMRecord::setAttribute(std::string &attributeName, const std::string& attributeValue) {
    ReadUtils::assertAttributeNameIsLegal(attributeName);
    setAttribute(attributeName, new std::string(attributeValue), String_Type, 0);
}

void SAMRecord::setAttribute(short tag, void *value, Void_Type type, int length) {
    setAttribute(tag, value, type, length, false);
}

void SAMRecord::setAttribute(std::string &tag, void *value, Void_Type type, int length) {
    setAttribute(SAMUtils::makeBinaryTag(tag), value, type, length);
}

int SAMRecord::getSoftStart() {
    std::shared_ptr<SAMRecord> ptr(new SAMRecord(*this));
    return ReadUtils::getSoftStart(ptr);
}

int SAMRecord::getSoftEnd() {
    std::shared_ptr<SAMRecord> ptr(new SAMRecord(*this));
    return ReadUtils::getSoftStart(ptr);
}

std::string &SAMRecord::getContig() {
    return getReadUnmappedFlag() ? (std::string&)"" : mReferenceName;
}

std::string SAMRecord::getAttributeAsString(std::string &attributeName) {
    ReadUtils::assertAttributeNameIsLegal(attributeName);
    if(mAttributes == nullptr) {
        return "";
    } else {
        SAMBinaryTagAndValue* tmp = mAttributes->find(SAMUtils::makeBinaryTag(attributeName));
        std::string ret;
        switch (tmp->type) {
            case Uint8_t_Array_Type:
            {
                char* val = (char*) tmp->value;
                if(tmp->length == 0) {
                    ret = "";
                    break;
                }
                char * newVal = new char[tmp->length+1];
                memcpy(newVal, val, tmp->length);
                newVal[tmp->length] = 0;
                ret = std::string(newVal);
                delete[] newVal;
                break;
            }
            case Uint8_Type:
            {
                char* val = (char*) tmp->value;
                ret = *val;
                break;
            }
            case Int_Array_Type:
            {
                std::stringstream ss;
                for(int i = 0; i < tmp->length; i++) {
                    ss << ((int*)tmp->value)[i];
                    ss >> ret;
                }
            }
            case Integer_Type:
            {
                std::stringstream ss;
                ss << *((int*)tmp->value);
                ss >> ret;
            }
            case Float_Array_Type:
            {
                std::stringstream ss;
                for(int i = 0; i < tmp->length; i++) {
                    ss << ((float *)tmp->value)[i];
                    ss >> ret;
                }
            }
            case Float_Type:
            {
                std::stringstream ss;
                ss << *((float *)tmp->value);
                ss >> ret;
            }
            case Long_Type:
            {
                std::stringstream ss;
                ss << *((long *)tmp->value);
                ss >> ret;
            }
            case String_Type:
            {
                return *((std::string*)tmp->value);
            }
            case Short_Type:
            {
                std::stringstream ss;
                ss << *((short *)tmp->value);
                ss >> ret;
            }
            case Short_Array_Type:
            {
                std::stringstream ss;
                for(int i = 0; i < tmp->length; i++) {
                    ss << ((short *)tmp->value)[i];
                    ss >> ret;
                }
            }
        }
        return ret;
    }
}

bool SAMRecord::isReverseStrand() const {
    return (mFlags & 16) != 0;
}

bool SAMRecord::mateIsReverseStrand() {
    Mutect2Utils::validateArg(isPaired(), "Cannot get mate information for an unpaired read");
    requireReadPaired();
    return getMateNegativeStrandFlagUnchecked();
}

bool SAMRecord::getMateNegativeStrandFlagUnchecked() {
    return (mFlags & 32) != 0;
}

bool SAMRecord::isFirstOfPair() {
    return isPaired() && getFirstOfPairFlag();
}

bool SAMRecord::getFirstOfPairFlag() {
    requireReadPaired();
    return (mFlags & 64) != 0;
}

bool SAMRecord::isSecondOfPair() {
    return isPaired() && getSecondOfPairFlag();
}

bool SAMRecord::getSecondOfPairFlag() {
    requireReadPaired();
    return (mFlags & 128) != 0;
}

bool SAMRecord::isSecondaryAlignment() const {
    return (mFlags & 256) != 0;
}

bool SAMRecord::failsVendorQualityCheck() const {
    return (mFlags & 512) != 0;
}

bool SAMRecord::isDuplicate() const {
    return (mFlags & 1024) != 0;
}

bool SAMRecord::isSupplementaryAlignment() const {
    return (mFlags & 2048) != 0;
}

SAMRecord::SAMRecord(bam1_t *read, SAMFileHeader* samFileHeader, bool load) {
    uint32_t * res = bam_get_cigar(read);
    uint32_t n = read->core.n_cigar;
    std::vector<CigarElement> nCigarElements;
    for(int i = 0; i < n; i++) {
        int length = (int) (res[i] >> 4);
        CigarOperator tmp_cigarOperator = CigarOperatorUtils::binaryToEnum((int) (res[i] & 0xf));
        nCigarElements.emplace_back(CigarElement(length, tmp_cigarOperator));
    }
    mCigar = new Cigar(nCigarElements);
    mFlags = read->core.flag;
    mMappingQuality = read->core.qual;
    mAlignmentStart = read->core.pos;
    mAlignmentEnd = mAlignmentStart + static_cast<int>(bam_cigar2rlen(n, res)) - 1;
    mReferenceName = std::string(samFileHeader->getSequenceDictionary().getSequences()[read->core.tid].getSequenceName());
    mMateReferenceName = std::string(samFileHeader->getSequenceDictionary().getSequences()[read->core.mtid].getSequenceName());
    mReadName = std::string(bam_get_qname(read));
    baseLength = read->core.l_qseq;
    baseQualitiesLength = static_cast<int>(bam_cigar2qlen(n, res));
    if(load) {
        uint8_t * bases = bam_get_seq(read);
        mReadBases = new uint8_t[baseLength+1]{0};
        mBaseQualities = new uint8_t[baseLength+1]{0};
        for(int i = 0; i < baseLength; i++) {
            uint8_t tmp = bam_seqi(bases, i);
            switch(tmp) {
                case 1 : {
                    mReadBases[i] = 'A' ;
                    break;
                }
                case 2 : {
                    mReadBases[i] = 'C';
                    break;
                }
                case 4 : {
                    mReadBases[i] = 'G';
                    break;
                }
                case 8 : {
                    mReadBases[i] = 'T';
                    break;
                }
                default : {
                    mReadBases[i] = 'N';
                    break;
                }
            }
        }
        uint8_t * qual = bam_get_qual(read);
        memcpy(mBaseQualities, qual, baseQualitiesLength);
    } else {
        mReadBases = nullptr;
        mBaseQualities = nullptr;
    }
    mMateAlignmentStart = static_cast<int>(read->core.mpos);
    mInferredInsertSize = static_cast<int>(read->core.isize);
    mAttributes = nullptr;
}

SAMRecord::~SAMRecord() {
    delete mCigar;
    delete mAttributes;
    //delete[] mReadBases;
    //delete[] mBaseQualities;
    //std::cout << "finished" << std::endl;
}

int SAMRecord::getAdaptorBoundary() {
    if(isCalAdaptorBoundary)
        return adaptorBoundary;
    std::shared_ptr<SAMRecord> ptr(new SAMRecord(*this));
    adaptorBoundary = ReadUtils::getAdaptorBoundary(ptr);
    isCalAdaptorBoundary = true;
    return adaptorBoundary;
}

SAMRecord::SAMRecord(const SAMRecord &other) : mFlags(other.mFlags), baseLength(other.baseLength), baseQualitiesLength(other.baseQualitiesLength),
mAlignmentStart(other.mAlignmentStart), mAlignmentEnd(other.mAlignmentEnd), mMateAlignmentStart(other.mMateAlignmentStart), mMappingQuality(other.mMappingQuality), mInferredInsertSize(other.mInferredInsertSize),
mReferenceName(other.mReferenceName), mMateReferenceName(other.mMateReferenceName), mReadName(other.mReadName){
    mAttributes = nullptr;
    if(other.mReadBases != nullptr){
        mReadBases = new uint8_t[baseLength+1]{0};
        memcpy(mReadBases, other.mReadBases, baseLength);
    }
    else{
        mReadBases = nullptr;
    }
    if(other.mBaseQualities != nullptr) {
        mBaseQualities = new uint8_t[baseQualitiesLength+1]{0};
        memcpy(mBaseQualities, other.mBaseQualities, baseQualitiesLength);
    }
    else{
        mBaseQualities = nullptr;
    }
    mCigar = new Cigar(other.mCigar->getCigarElements());
}

SimpleInterval SAMRecord::getLoc() {
    return {mReferenceName, mAlignmentStart, mAlignmentEnd};
}

int SAMRecord::getEndAfterFliter() {
    return mAlignmentEnd;
}



