//
// Created by lhh on 10/19/21.
//

#include "ReadFilter.h"
#include "Mutect2Engine.h"
#include "read/CigarUtils.h"
#include "QualityUtils.h"
#include "read/ReadUtils.h"

bool ReadFilter::ReadLengthTest()
{
    if(originalRead.getLength() >= Mutect2Engine::MIN_READ_LENGTH)
        return true;
    else
        return false;
}

bool ReadFilter::NotSecondaryAlignmentTest() const {
    return ! originalRead.isSecondaryAlignment();
}

bool ReadFilter::GoodCigarTest() {
    return CigarUtils::isGood(originalRead.getCigar());
}

ReadFilter::ReadFilter(bam1_t *read, SAMFileHeader* header) : originalRead(read, header, false), header(header){
}

bool ReadFilter::NonZeroReferenceLengthAlignmentTest() {
    for(CigarElement element : originalRead.getCigarElements()) {
        if(CigarOperatorUtils::getConsumesReferenceBases(element.getOperator()) && element.getLength() > 0) {
            return true;
        }
    }
    return false;
}

bool ReadFilter::PassesVendorQualityCheck() const {
    return ! originalRead.failsVendorQualityCheck();
}

bool ReadFilter::MappedTest() {
    return ! originalRead.isUnmapped();
}

bool ReadFilter::MappingQualityAvailableTest() {
    return originalRead.getMappingQuality() != QualityUtils::MAPPING_QUALITY_UNAVALIABLE;
}

bool ReadFilter::NotDuplicateTest() {
    return ! originalRead.isDuplicate();
}

bool ReadFilter::MappingQualityTest() {
    return originalRead.getMappingQuality() >= 20;
}

bool ReadFilter::MappingQualityNotZeroTest() {
    return originalRead.getMappingQuality() != 20;
}

bool ReadFilter::WellformedTest() {
    return (originalRead.isUnmapped() || originalRead.getStart() > 0) &&
            (originalRead.isUnmapped() || (originalRead.getEnd() - originalRead.getStart() + 1) >= 0) &&
            ReadUtils::alignmentAgreesWithHeader(header, &originalRead) &&
            // ! originalRead.getReadGroup().empty() &&
            originalRead.getLength() == originalRead.getBaseQualitiesLength() &&
            (originalRead.isUnmapped() || originalRead.getLength() == Cigar::getReadLength(originalRead.getCigarElements())) &&
            (originalRead.getLength() > 0) &&
            (! CigarUtils::containsNOperator(originalRead.getCigarElements()));
}

bool ReadFilter::test() {
    bool ret1 = ReadLengthTest();
    bool ret2 = NonZeroReferenceLengthAlignmentTest();
    bool ret3 = NotDuplicateTest();
    bool ret4 = NotSecondaryAlignmentTest();
    bool ret5 = GoodCigarTest();
    bool ret6 = PassesVendorQualityCheck();
    bool ret7 = MappedTest();
    bool ret8 = MappingQualityAvailableTest();
    bool ret9 = MappingQualityNotZeroTest();
    bool ret10 = MappingQualityTest();
    bool ret11 = WellformedTest();
    return ReadLengthTest()&& NonZeroReferenceLengthAlignmentTest() && NotDuplicateTest() && NotSecondaryAlignmentTest() &&
    GoodCigarTest() && PassesVendorQualityCheck() && MappedTest() && MappingQualityAvailableTest() && MappingQualityNotZeroTest() && MappingQualityTest() && WellformedTest();

//    return true;
}
