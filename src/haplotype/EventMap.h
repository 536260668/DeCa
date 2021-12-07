//
// Created by 梦想家xixi on 2021/12/4.
//

#ifndef MUTECT2CPP_MASTER_EVENTMAP_H
#define MUTECT2CPP_MASTER_EVENTMAP_H

#include "Allele.h"
#include "Haplotype.h"
#include "variantcontext/VariantContext.h"

class EventMap {
private:
    static const int MAX_EVENTS_PER_HAPLOTYPE = 3;
    static const int MAX_INDELS_PER_HAPLOTYPE = 2;
    Haplotype* haplotype;
    uint8_t * ref;
    int refLength;
    Locatable* refLoc;
    std::string sourceNameToAdd;
    std::map<int, VariantContext*> variantMap;

public:
    static const Allele* SYMBOLIC_UNASSEMBLED_EVENT_ALLELE;
    EventMap(Haplotype* haplotype, uint8_t* ref, int refLength, Locatable* refLoc, std::string sourceNameToAdd, int maxMnpDistance);
    void addVC(VariantContext* vc, bool merge);

protected:
    static const int MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION = 3;
    void processCigarForInitialEvents(int maxMnpDistance);
    VariantContext* makeBlock(VariantContext* vc1, VariantContext* vc2);
};


#endif //MUTECT2CPP_MASTER_EVENTMAP_H
