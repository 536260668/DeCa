//
// Created by 梦想家xixi on 2021/12/18.
//

#include "ReadUtils.h"

SAMRecord *ReadUtils::emptyRead(SAMRecord *read) {
    SAMRecord* emptyRead = new SAMRecord(*read);
    emptyRead->setIsUnmapped();
    emptyRead->setMappingQuality(0);
    emptyRead->setCigar(new Cigar());
    emptyRead->setBases(nullptr, 0);
    emptyRead->setBaseQualities(nullptr, 0);
    emptyRead->clearAttributes();

    std::string readGroup = read->getReadGroup();
    if(readGroup.empty()) {
        emptyRead->setAttribute((std::string&)"RG", readGroup);
    }
    return emptyRead;
}

void ReadUtils::assertAttributeNameIsLegal(std::string &attributeName) {
    if(attributeName.empty() || attributeName.length() !=2) {
        throw std::invalid_argument("Read attribute invalid: attribute names must be non-null two-character Strings matching the pattern /[A-Za-z][A-Za-z0-9]/");
    }
}
