//
// Created by 梦想家xixi on 2022/1/15.
//

#include "ReferenceConfidenceModel.h"

#include <utility>
#include "Mutect2Utils.h"

std::shared_ptr<Haplotype>
ReferenceConfidenceModel::createReferenceHaplotype(std::shared_ptr<AssemblyRegion> activeRegion, std::shared_ptr<uint8_t[]> refBases, int & length,
                                                   SimpleInterval &paddedReferenceLoc) {
   int alignmentStart = activeRegion->getExtendedSpan().getStart() - paddedReferenceLoc.getStart();
   if(alignmentStart < 0) {
       throw std::invalid_argument("Bad alignment start in createReferenceHaplotype");
   }
    std::shared_ptr<Haplotype> refHaplotype(new Haplotype(std::move(refBases), length, true));
    refHaplotype->setAlignmentStartHapwrtRef(alignmentStart);
    std::shared_ptr<Cigar> c(new Cigar());
    c->add(CigarElement(refHaplotype->getBasesLength(), M));
    refHaplotype->setCigar(c);
    return refHaplotype;
}
