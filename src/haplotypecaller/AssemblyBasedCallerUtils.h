//
// Created by 梦想家xixi on 2021/12/8.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYBASEDCALLERUTILS_H
#define MUTECT2CPP_MASTER_ASSEMBLYBASEDCALLERUTILS_H

#include "haplotype/Haplotype.h"
#include "AssemblyRegion.h"
#include "ReferenceCache.h"
#include "Mutect2/AssemblyResultSet.h"
#include "M2ArgumentCollection.h"
#include "ReadThreadingAssembler.h"

class AssemblyBasedCallerUtils {
public:
    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param region the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the interval which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    static std::shared_ptr<Haplotype> createReferenceHaplotype(AssemblyRegion & region, SimpleInterval referencePadding, ReferenceCache & cache);

    static std::shared_ptr<AssemblyResultSet> assembleReads(AssemblyRegion& region, M2ArgumentCollection & argumentCollection, SAMFileHeader* header, ReferenceCache & cache, ReadThreadingAssembler& assemblyEngine);

    static const int REFERENCE_PADDING_FOR_ASSEMBLY = 500;

    static SimpleInterval getPaddedReferenceLoc(AssemblyRegion& region, int referencePadding, SAMFileHeader* header);
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYBASEDCALLERUTILS_H
