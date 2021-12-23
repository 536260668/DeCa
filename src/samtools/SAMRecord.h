//
// Created by 梦想家xixi on 2021/12/20.
//

#ifndef MUTECT2CPP_MASTER_SAMRECORD_H
#define MUTECT2CPP_MASTER_SAMRECORD_H

#include <string>
#include "Locatable.h"
#include "Cigar.h"
#include "SAMBinaryTagAndValue.h"

class SAMRecord{
private:
    uint8_t * mReadBases;
    int baseLength;
    uint8_t * mBaseQualities;
    int baseQualitiesLength;
    std::string mReadName;
    std::string mReferenceName;
    int mAlignmentStart;
    int mAlignmentEnd;
    int mMappingQuality;
    Cigar* mCigar;
    int mFlags;
    std::string mMateReferenceName;
    int mMateAlignmentStart;
    int mInferredInsertSize;
    SAMBinaryTagAndValue* mAttributes;
    void setFlag(bool flag, int bit);
    void requireReadPaired();
    bool getMateUnmappedFlagUnchecked();

public:
    SAMRecord(uint8_t* base, int baseLength, uint8_t* baseQualities, int baseQualitiesLength, std::string &name);
    static const std::string NO_ALIGNMENT_REFERENCE_NAME;
    static const int NO_ALIGNMENT_START = 0;
    static const int NO_MAPPING_QUALITY = 0;
    std::string& getName();
    void setName(std::string& name);
    int getLength() const;
    bool isEmpty() const {return getLength() == 0;}
    void setPosition(std::string& contig, int start);
    void setPosition(Locatable* locatable);
    std::string & getAssignedContig();
    int getAssignedStart();
    int getUnclippedStart();
    int getUnclippedEnd();
    bool getReadUnmappedFlag();
    bool isUnmapped();
    bool mateIsUnmapped();
    bool getMateUnmappedFlag();
    bool isPaired();
    std::string& getMateContig();
    int getMateStart();
    void setMatePosition(std::string& contig, int start);
    void setIsPaired(bool isPaired);
    void setReadPairedFlag(bool flag);
    void setProperPairFlag(bool flag);
    void setMatePosition(Locatable* locatable);
    int getFragmentLength() const;
    void setFragmentLength(int fragmentLength);
    int getMappingQuality() const;
    void setMappingQuality(int mappingQuality);
    uint8_t * getBases();
    uint8_t * getBasesNoCopy();
    uint8_t getBase(const int i) {return getBasesNoCopy()[i];}
    int getLength();
    void setBases(uint8_t *bases, int length);
    uint8_t * getBaseQualities();
    uint8_t * getBaseQualitiesNoCopy();
    int getBaseQualitiesLength();
    uint8_t getBaseQuality(const int i) {return getBaseQualities()[i];}
    void setBaseQualities(uint8_t* baseQualities, int length);
    Cigar* getCigar();
    std::vector<CigarElement>& getCigarElements();
    CigarElement getCigarElement(int index);
    void setCigar(Cigar* cigar);
    bool getProperPairFlag();
    bool getProperPairFlagUnchecked() const;
    bool isProperlyPaired();
    void setIsProperlyPaired(bool isProperlyPaired);
    int getStart();
    int getEnd();
    void setIsUnmapped();
    void setReadUnmappedFlag(bool flag);
    void clearAttributes();
    void* getAttribute(short tag);
    std::string& getReadGroup();
    void setAttribute(std::string& attributeName, const std::string& attributeValue);
    void setAttribute(std::string &tag, void* value, Void_Type type, int length);
    void setAttribute(short tag, void* value, Void_Type type, int length);
    int getSoftStart();
    int getSoftEnd();
    std::string & getContig();
    std::string getAttributeAsString(std::string & attributeName);
    bool isReverseStrand() const;
    bool mateIsReverseStrand();
    bool getMateNegativeStrandFlagUnchecked();
    bool isFirstOfPair();
    bool getFirstOfPairFlag();
    bool isSecondOfPair();
    bool getSecondOfPairFlag();
    bool isSecondaryAlignment() const;
    bool failsVendorQualityCheck() const;
    bool isDuplicate() const;
    bool isSupplementaryAlignment() const;

private:
    void setAttribute(short tag, void* value, Void_Type type, int length, bool isUnsignedArray);
};


#endif //MUTECT2CPP_MASTER_SAMRECORD_H
