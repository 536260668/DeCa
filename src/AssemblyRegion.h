//
// Created by 梦想家xixi on 2021/10/14.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYREGION_H
#define MUTECT2CPP_MASTER_ASSEMBLYREGION_H

#include <vector>
#include "SimpleInterval.h"
#include "ActivityProfileState.h"

class AssemblyRegion {
private:
    //SAMFileHeader header;

    /**
     * The reads included in this assembly region.  May be empty upon creation, and expand / contract
     * as reads are added or removed from this region.
     */
    //std::vector<GATKReads> reads;

    /**
     * An ordered list (by genomic coordinate) of the ActivityProfileStates that went
     * into this assembly region.  May be empty, which says that no supporting states were
     * provided when this region was created.
     */
     std::vector<ActivityProfileState> supportingStates;

    /**
    * The raw span of this assembly region, not including the region extension
    */
    SimpleInterval activeRegionLoc;

    /**
     * The span of this assembly region on the genome, including the region extension
     */
     SimpleInterval extendedLoc;

    /**
    * The extension, in bp, of this region. The extension is >= 0 bp in size, and indicates how much padding was
    * requested for the region.
    */
    const int extension;

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
     bool isActivate;

    /**
    * The span of this assembly region, including the bp covered by all reads in this
    * region.  This union of extensionLoc and the loc of all reads in this region.
    *
    * Must be at least as large as extendedLoc, but may be larger when reads
    * partially overlap this region.
    */
    SimpleInterval spanIncludingReads;

    /**
     * Indicates whether the region has been finalized
     */
     bool hasBeenFinalized;

    SimpleInterval* trimIntervalToContig(std::string& contig, int start, int stop);

    void checkStates(SimpleInterval& activeRegion);

public:
    AssemblyRegion(SimpleInterval& activeRegionLoc, std::vector<ActivityProfileState> supportingStates, bool isActivate, int extension);
    virtual ~AssemblyRegion();
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYREGION_H
