//
// Created by 梦想家xixi on 2021/12/8.
//

#include "AssemblyBasedCallerUtils.h"
#include "haplotypecaller/ReferenceConfidenceModel.h"
#include "clipping/ReadClipper.h"
#include "read/ReadUtils.h"

std::shared_ptr<Haplotype>
AssemblyBasedCallerUtils::createReferenceHaplotype(std::shared_ptr<AssemblyRegion> region, SimpleInterval &referencePadding,
                                                   ReferenceCache &cache) {
    int length = 0;
    std::shared_ptr<uint8_t[]> tmp = region->getAssemblyRegionReference(&cache, 0, length);
    return ReferenceConfidenceModel::createReferenceHaplotype(region, tmp, length, referencePadding);
}

std::shared_ptr<AssemblyResultSet>
AssemblyBasedCallerUtils::assembleReads(std::shared_ptr<AssemblyRegion> region, M2ArgumentCollection &argumentCollection,
                                        SAMFileHeader *header, ReferenceCache &cache,
                                        ReadThreadingAssembler &assemblyEngine) {
    finalizeRegion(region, false, false, 9, header, false);
    int refLength = 0;
    std::shared_ptr<uint8_t[]> fullReferenceWithPadding = region->getAssemblyRegionReference(&cache, REFERENCE_PADDING_FOR_ASSEMBLY, refLength);
    SimpleInterval paddedReferenceLoc = getPaddedReferenceLoc(region, REFERENCE_PADDING_FOR_ASSEMBLY, header);
    std::shared_ptr<Haplotype> refHaplotype = createReferenceHaplotype(region, paddedReferenceLoc, cache);
    std::shared_ptr<AssemblyResultSet> assemblyResultSet = assemblyEngine.runLocalAssembly(region, refHaplotype, fullReferenceWithPadding, refLength, & paddedReferenceLoc,
                                                                          nullptr);
    return assemblyResultSet;
}

SimpleInterval
AssemblyBasedCallerUtils::getPaddedReferenceLoc(std::shared_ptr<AssemblyRegion> region, int referencePadding, SAMFileHeader *header) {
    int padLeft = std::max(region->getExtendedSpan().getStart() - referencePadding, 1);
    int padRight = std::min(region->getExtendedSpan().getEnd() + referencePadding, header->getSequenceDictionary().getSequence(region->getExtendedSpan().getContig()).getSequenceLength());
    return {region->getExtendedSpan().getContig(), padLeft, padRight};
}

void
AssemblyBasedCallerUtils::finalizeRegion(const std::shared_ptr<AssemblyRegion>& region, bool errorCorrectReads, bool dontUseSoftClippedBases,
                                         uint8_t minTailQuality, SAMFileHeader *header,
                                         bool correctOverlappingBaseQualities) {
    if(region->isFinalized())
        return;
    std::vector<std::shared_ptr<SAMRecord>> readsToUse;
    for(const std::shared_ptr<SAMRecord>& myRead : region->getReads()) {
        uint8_t minTailQualityToUse = errorCorrectReads ? 6 : minTailQuality;
        std::shared_ptr<SAMRecord> clippedRead = ReadClipper::hardClipLowQualEnds(myRead, minTailQuality);
        clippedRead = dontUseSoftClippedBases || ! ReadUtils::hasWellDefinedFragmentSize(clippedRead) ?
                ReadClipper::hardClipSoftClippedBases(clippedRead) : ReadClipper::revertSoftClippedBases(clippedRead);

        clippedRead = clippedRead->isUnmapped() ? clippedRead : ReadClipper::hardClipAdaptorSequence(clippedRead);
        if(! clippedRead->isEmpty() && clippedRead->getCigar()->getReadLength() > 0 ) {
            clippedRead = ReadClipper::hardClipToRegion( clippedRead, region->getExtendedSpan().getStart(), region->getExtendedSpan().getEnd());
            if ( region->readOverlapsRegion(clippedRead) && clippedRead->getLength() > 0 ) {
                readsToUse.emplace_back((clippedRead == myRead) ? std::shared_ptr<SAMRecord>(new SAMRecord(*clippedRead)) : clippedRead);
            }
        }
    }
    region->clearReads();
    region->addAll(readsToUse);
    region->setFinalized(true);
}
