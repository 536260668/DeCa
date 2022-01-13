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
    static bool ReadLengthTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool NotSecondaryAlignmentTest(std::shared_ptr<SAMRecord> & originalRead) ;
    static bool GoodCigarTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool NonZeroReferenceLengthAlignmentTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool PassesVendorQualityCheck(std::shared_ptr<SAMRecord> &originalRead);
    static bool MappedTest(std::shared_ptr<SAMRecord> & originalRead );
    static bool MappingQualityAvailableTest(std::shared_ptr<SAMRecord> &originalRead);
    static bool NotDuplicateTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool MappingQualityTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool MappingQualityNotZeroTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool WellformedTest(std::shared_ptr<SAMRecord> & originalRead, SAMFileHeader* header);
    static bool test(std::shared_ptr<SAMRecord> & originalRead, SAMFileHeader* header);

private:
};


#endif //MUTECT2CPP_MASTER_READFILTER_H
