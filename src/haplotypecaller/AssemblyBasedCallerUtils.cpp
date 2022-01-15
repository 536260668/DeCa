//
// Created by 梦想家xixi on 2021/12/8.
//

#include "AssemblyBasedCallerUtils.h"
#include "haplotypecaller/ReferenceConfidenceModel.h"

std::shared_ptr<Haplotype>
AssemblyBasedCallerUtils::createReferenceHaplotype(AssemblyRegion &region, SimpleInterval referencePadding,
                                                   ReferenceCache &cache) {
    return ReferenceConfidenceModel::createReferenceHaplotype(region, region.getAssemblyRegionReference(&cache, 0), referencePadding);
}
