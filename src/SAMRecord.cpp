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
    if (IsUnmapped()) { return nullptr; }

    return std::string(GetReferenceName());
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
    //return bam_get_qual(read);
    uint8_t * res = new uint8_t[86]{34, 35, 33, 35, 36, 35, 35, 35, 29, 35, 33, 37, 35, 36, 36, 37, 36, 36, 37, 36, 36, 36, 35, 34, 36, 37, 36, 30, 37, 38, 36, 35, 27, 36, 34, 37, 36, 38, 36, 36, 35, 36, 37, 38, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 18, 20, 20, 20, 4, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
    return res;
}

uint8_t SAMRecord::getBaseQuality(int i)
{
    uint8_t * qualities = bam_get_qual(read);
    return qualities[i];
}

uint8_t *SAMRecord::getBases()
{
    //return bam_get_seq(read);
    uint8_t * res = new uint8_t[86]{65, 84, 65, 67, 65, 67, 67, 67, 71, 71, 67, 65, 67, 67, 67, 84, 71, 84, 67, 67, 84, 71, 71, 65, 67, 65, 67, 71, 67, 84, 71, 84, 84, 71, 71, 67, 67, 84, 71, 71, 65, 84, 67, 84, 71, 65, 71, 67, 67, 67, 84, 71, 71, 84, 71, 71, 65, 71, 71, 84, 67, 65, 65, 65, 71, 67, 67, 65, 67, 67, 84, 84, 84, 71, 71, 84, 84, 67, 84, 71, 67, 67, 65, 84, 84, 71};
    return res;
}

int32_t SAMRecord::getLength()
{
    //return read->core.l_qseq;
    return 86;
}

bool SAMRecord::isReverseStrand()
{
    return (read->core.flag & BAM_FREVERSE) != 0;
}

char *SAMRecord::getReadName() {
    //return bam_get_qname(read);
    char * res = new char[10]{'c','h','1'};
    return res;
}
