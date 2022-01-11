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
    static bool ReadLengthTest(SAMRecord & originalRead);
    static bool NotSecondaryAlignmentTest(SAMRecord & originalRead) ;
    static bool GoodCigarTest(SAMRecord & originalRead);
    static bool NonZeroReferenceLengthAlignmentTest(SAMRecord & originalRead);
    static bool PassesVendorQualityCheck(SAMRecord & originalRead);
    static bool MappedTest(SAMRecord & originalRead );
    static bool MappingQualityAvailableTest(SAMRecord & originalRead);
    static bool NotDuplicateTest(SAMRecord & originalRead);
    static bool MappingQualityTest(SAMRecord & originalRead);
    static bool MappingQualityNotZeroTest(SAMRecord & originalRead);
    static bool WellformedTest(SAMRecord & originalRead, SAMFileHeader* header);
    static bool test(SAMRecord & originalRead, SAMFileHeader* header);

private:
};


#endif //MUTECT2CPP_MASTER_READFILTER_H
