//
// The class used to filter out reads for Mutect2Cpp-master
// Created by lhh on 10/19/21.
//

#ifndef MUTECT2CPP_MASTER_READFILTER_H
#define MUTECT2CPP_MASTER_READFILTER_H

#include "htslib/sam.h"
#include "samtools/SAMRecord.h"

class ReadFilter {  // TODO: add the filter left
public:
    bool ReadLengthTest();
    bool NotSecondaryAlignmentTest() const;
    bool GoodCigarTest();
    bool NonZeroReferenceLengthAlignmentTest();
    bool PassesVendorQualityCheck() const;
    bool MappedTest();
    bool MappingQualityAvailableTest();
    bool NotDuplicateTest();
    bool MappingQualityTest();
    bool MappingQualityNotZeroTest();
    bool WellformedTest();
    bool test();
    explicit ReadFilter(bam1_t* read, SAMFileHeader* header);

private:
    SAMRecord originalRead;
    SAMFileHeader* header;
};


#endif //MUTECT2CPP_MASTER_READFILTER_H
