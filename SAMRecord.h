//
// The class that represent a sam record based on htslib
// Created by lhh on 10/14/21.
//

#ifndef MUTECT2CPP_MASTER_SAMRECORD_H
#define MUTECT2CPP_MASTER_SAMRECORD_H

#include <string>
#include "htslib/sam.h"

#define MAPPING_QUALITY_UNAVAILABLE 255
#define NO_ALIGNMENT_REFERENCE_NAME "*"
#define NO_ALIGNMENT_START 0
#define UNSET_POSITION 0

class SAMRecord {
private:
    bam1_t * read;
    sam_hdr_t * header;

    //test
    uint8_t * bases;
    int baseLength;
    uint8_t * baseQualities;
    int baseQualitiesLength;
    std::string name;

public:
    SAMRecord(bam1_t * samRecord, sam_hdr_t * header);
    SAMRecord(uint8_t* base, int baseLength, uint8_t* baseQualities, int baseQualitiesLength, std::string &name);
    SAMRecord(const SAMRecord & other);
    /**
    * zero-based start
    */
    hts_pos_t GetStart();

    /**
     * Does the read have a position assigned to it for sorting purposes.
     * @return `true if this read has no assigned position or contig.
     */
    bool IsUnmapped();

    /**
     * Only to encapsulate sam_hdr_tid2name function in sam.h
     */
    const char* GetReferenceName();

    std::string getContig();

    /**
     * @return zero-based end
     */
    hts_pos_t GetEnd();

    std::string getMateContig();

    const char*  getMateReferenceName();

    hts_pos_t getFragmentLength();

    uint8_t * getBaseQualities();

    uint8_t getBaseQuality(int i);

    /**
     * @return pointer to sequence
     * Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
     * 8 for T and 15 for N. Two bases are packed in one byte with the base
     * at the higher 4 bits having smaller coordinate on the read.
     */
     uint8_t * getBases();

     int32_t getLength();

     /**
     * @return True if this read is on the reverse strand as opposed to the forward strand, otherwise false.
     */
     bool isReverseStrand();

     std::string getReadName();
};


#endif //MUTECT2CPP_MASTER_SAMRECORD_H
