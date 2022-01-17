//
// Created by lhh on 10/14/21.
//

#include "SAMRecord.h"

SAMRecord::SAMRecord(bam1_t *samRecord, sam_hdr_t *header):read(samRecord), header(header)
{}

hts_pos_t SAMRecord::GetStart()
{
    return read->core.pos;
}

bool SAMRecord::IsUnmapped()
{
    const char * refName = GetReferenceName();
    return (read->core.flag & BAM_FUNMAP) != 0 ||
    strcmp(refName, NO_ALIGNMENT_REFERENCE_NAME) == 0 ||
    read->core.pos == NO_ALIGNMENT_START;
}

const char *SAMRecord::GetReferenceName()
{
    return sam_hdr_tid2name(header, read->core.tid);
}

std::string SAMRecord::getContig()
{
    //if (IsUnmapped()) { return nullptr; }
    return "chr1";
    //return std::string(GetReferenceName());
}

hts_pos_t SAMRecord::GetEnd()
{
    return bam_endpos(read) - 1;
}

std::string SAMRecord::getMateContig()
{
    return std::string(getMateReferenceName());
}

const char* SAMRecord::getMateReferenceName()
{
    return sam_hdr_tid2name(header, read->core.mtid);
}

hts_pos_t SAMRecord::getFragmentLength()
{
    return read->core.isize;
}

uint8_t * SAMRecord::getBaseQualities()
{
    return baseQualities;
    //return bam_get_qual(read);
}

uint8_t SAMRecord::getBaseQuality(int i)
{
    uint8_t * qualities = bam_get_qual(read);
    return qualities[i];
}

uint8_t *SAMRecord::getBases()
{
    return bases;
    //return bam_get_seq(read);
}

int32_t SAMRecord::getLength()
{
    return baseLength;
    // return read->core.l_qseq;
}

bool SAMRecord::isReverseStrand()
{
    return (read->core.flag & BAM_FREVERSE) != 0;
}

std::string SAMRecord::getReadName() {
    return name.c_str();
    //return bam_get_qname(read);
}

SAMRecord::SAMRecord(uint8_t *base, int baseLength, uint8_t *baseQualities, int baseQualitiesLength,
                     std::string &name) : bases(base), baseLength(baseLength), baseQualities(baseQualities),baseQualitiesLength(baseQualitiesLength), name(name), header(
        nullptr), read(nullptr){}

SAMRecord::SAMRecord(const SAMRecord &other) {
    read = other.read;
    header = other.header;
    bases = other.bases;
    baseLength = other.baseLength;
    baseQualities = other.baseQualities;
    baseQualitiesLength = other.baseQualitiesLength;
}
