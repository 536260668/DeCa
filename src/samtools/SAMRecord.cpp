//
// Created by 梦想家xixi on 2021/12/20.
//

#include "SAMRecord.h"
#include "SAMFlag.h"
#include "ReadConstants.h"
#include "SAMUtils.h"
#include "ReadUtils.h"

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
    if(!contig.empty() || contig == NO_ALIGNMENT_REFERENCE_NAME || start < 1) {
        throw std::invalid_argument("contig must be non-null and start must be >= 1");
    }
    mReferenceName = contig;
    mAlignmentStart = start;
    mAlignmentEnd = 0;
    setFlag(false, SAMFlag::READ_UNMAPPED().flag);
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
    return (mFlags & SAMFlag::READ_UNMAPPED().flag) != 0;
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
    return (mFlags & SAMFlag::READ_PAIRED().flag) != 0;
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
    if((mFlags & SAMFlag::READ_PAIRED().flag) == 0) {
        throw std::invalid_argument("Inappropriate call if not paired read");
    }
}

bool SAMRecord::getMateUnmappedFlagUnchecked() {
    return (mFlags & SAMFlag::MATE_UNMAPPED().flag) != 0;
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
    setFlag(false, SAMFlag::MATE_UNMAPPED().flag);
}

void SAMRecord::setReadPairedFlag(bool flag) {
    setFlag(flag, SAMFlag::READ_PAIRED().flag);
}

void SAMRecord::setProperPairFlag(bool flag) {
    setFlag(flag, SAMFlag::PROPER_PAIR().flag);
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
    Mutect2Utils::validateArg(mappingQuality < 0 || mappingQuality > 255, "mapping quality must be >= 0 and <= 255");
    mMappingQuality = mappingQuality;
}

uint8_t *SAMRecord::getBases() {
    if(mReadBases != nullptr){
        uint8_t * ret = new uint8_t[baseLength];
        memcpy(ret, mReadBases, baseLength);
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
        uint8_t * ret = new uint8_t[baseQualitiesLength];
        memcpy(ret, mBaseQualities, baseQualitiesLength);
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
    mCigar = cigar;
    mAlignmentEnd = 0;
}

bool SAMRecord::getProperPairFlag() {
    requireReadPaired();
    return getProperPairFlagUnchecked();
}

bool SAMRecord::getProperPairFlagUnchecked() const {
    return (mFlags & SAMFlag::PROPER_PAIR().flag) != 0;
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
    setFlag(flag, SAMFlag::READ_UNMAPPED().flag);
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



