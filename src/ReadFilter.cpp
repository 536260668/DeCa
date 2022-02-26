//
// Created by lhh on 10/19/21.
//

#include "ReadFilter.h"
#include "Mutect2Engine.h"
#include "read/CigarUtils.h"
#include "QualityUtils.h"
#include "read/ReadUtils.h"


bool ReadFilter::NotSecondaryAlignmentTest(std::shared_ptr<SAMRecord> & originalRead) {
    return ! originalRead->isSecondaryAlignment();
}

bool ReadFilter::GoodCigarTest(std::shared_ptr<SAMRecord> & originalRead) {
    return CigarUtils::isGood(originalRead->getCigar());
}


bool ReadFilter::NonZeroReferenceLengthAlignmentTest(std::shared_ptr<SAMRecord> & originalRead) {
    for(const CigarElement& element : originalRead->getCigarElements()) {
        if(CigarOperatorUtils::getConsumesReferenceBases(element.getOperator()) && element.getLength() > 0) {
            return true;
        }
    }
    return false;
}

bool ReadFilter::PassesVendorQualityCheck(std::shared_ptr<SAMRecord> & originalRead) {
    return ! originalRead->failsVendorQualityCheck();
}

bool ReadFilter::MappedTest(std::shared_ptr<SAMRecord> & originalRead) {
    return ! originalRead->isUnmapped();
}

bool ReadFilter::MappingQualityAvailableTest(std::shared_ptr<SAMRecord> & originalRead) {
    return originalRead->getMappingQuality() != QualityUtils::MAPPING_QUALITY_UNAVALIABLE;
}

bool ReadFilter::NotDuplicateTest(std::shared_ptr<SAMRecord> & originalRead) {
    return ! originalRead->isDuplicate();
}

bool ReadFilter::MappingQualityTest(std::shared_ptr<SAMRecord> & originalRead) {
    return originalRead->getMappingQuality() >= 20;
}

bool ReadFilter::MappingQualityNotZeroTest(std::shared_ptr<SAMRecord> & originalRead) {
    return originalRead->getMappingQuality() != 0;
}

bool ReadFilter::WellformedTest(std::shared_ptr<SAMRecord> & originalRead, SAMFileHeader* header) {
    return (originalRead->isUnmapped() || originalRead->getStart() > 0) &&
            (originalRead->isUnmapped() || (originalRead->getEnd() - originalRead->getStart() + 1) >= 0) &&
            ReadUtils::alignmentAgreesWithHeader(header, originalRead) &&
            // ! originalRead.getReadGroup().empty() &&
            originalRead->getLength() == originalRead->getBaseQualitiesLength() &&
            (originalRead->isUnmapped() || originalRead->getLength() == Cigar::getReadLength(originalRead->getCigarElements())) &&
            (originalRead->getLength() > 0) &&
            (! CigarUtils::containsNOperator(originalRead->getCigarElements()));
}

bool ReadFilter::test(std::shared_ptr<SAMRecord> & originalRead, SAMFileHeader* header) {
//    if(originalRead.getStart() == 10056)
//        std::cout << "hello";
//    bool ret1 = ReadLengthTest();
//    bool ret2 = NonZeroReferenceLengthAlignmentTest();
//    bool ret3 = NotDuplicateTest();
//    bool ret4 = NotSecondaryAlignmentTest();
//    bool ret5 = GoodCigarTest();
//    bool ret6 = PassesVendorQualityCheck();
//    bool ret7 = MappedTest();
//    bool ret8 = MappingQualityAvailableTest();
//    bool ret9 = MappingQualityNotZeroTest();
//    bool ret10 = MappingQualityTest();
//    bool ret11 = WellformedTest();
    return ReadLengthTest(originalRead)&& NonZeroReferenceLengthAlignmentTest(originalRead) && NotDuplicateTest(originalRead) && NotSecondaryAlignmentTest(originalRead) &&
    GoodCigarTest(originalRead) && PassesVendorQualityCheck(originalRead) && MappedTest(originalRead) && MappingQualityAvailableTest(originalRead) && MappingQualityNotZeroTest(originalRead) && MappingQualityTest(originalRead) && WellformedTest(originalRead,header);

//    return true;
}

bool ReadFilter::ReadLengthTest(std::shared_ptr<SAMRecord> & originalRead) {
    return originalRead->getLength() > 30 && originalRead->getLength() < 2147483647;
}
