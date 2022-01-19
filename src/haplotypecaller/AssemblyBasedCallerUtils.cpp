//
// Created by 梦想家xixi on 2021/12/8.
//

#include "AssemblyBasedCallerUtils.h"
#include "haplotypecaller/ReferenceConfidenceModel.h"

std::shared_ptr<Haplotype>
AssemblyBasedCallerUtils::createReferenceHaplotype(AssemblyRegion &region, SimpleInterval &referencePadding,
                                                   ReferenceCache &cache) {
    int length = 0;
    uint8_t * tmp = region.getAssemblyRegionReference(&cache, 0, length);
    return ReferenceConfidenceModel::createReferenceHaplotype(region, tmp, length, referencePadding);
}

std::shared_ptr<AssemblyResultSet>
AssemblyBasedCallerUtils::assembleReads(AssemblyRegion &region, M2ArgumentCollection &argumentCollection,
                                        SAMFileHeader *header, ReferenceCache &cache,
                                        ReadThreadingAssembler &assemblyEngine) {
    int refLength = 0;
    uint8_t * fullReferenceWithPadding = region.getAssemblyRegionReference(&cache, REFERENCE_PADDING_FOR_ASSEMBLY, refLength);
    SimpleInterval paddedReferenceLoc = getPaddedReferenceLoc(region, REFERENCE_PADDING_FOR_ASSEMBLY, header);
    std::shared_ptr<Haplotype> refHaplotype = createReferenceHaplotype(region, paddedReferenceLoc, cache);
    std::shared_ptr<AssemblyResultSet> assemblyResultSet = assemblyEngine.runLocalAssembly(region, refHaplotype, fullReferenceWithPadding, refLength, & paddedReferenceLoc,
                                                                          nullptr);
    return assemblyResultSet;
}

SimpleInterval
AssemblyBasedCallerUtils::getPaddedReferenceLoc(AssemblyRegion &region, int referencePadding, SAMFileHeader *header) {
    int padLeft = std::max(region.getExtendedSpan().getStart() - referencePadding, 1);
    int padRight = std::min(region.getExtendedSpan().getEnd() + referencePadding, header->getSequenceDictionary().getSequence(region.getExtendedSpan().getContig()).getSequenceLength());
    return {region.getExtendedSpan().getContig(), padLeft, padRight};
}
