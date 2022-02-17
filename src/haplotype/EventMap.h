//
// Created by 梦想家xixi on 2021/12/4.
//

#ifndef MUTECT2CPP_MASTER_EVENTMAP_H
#define MUTECT2CPP_MASTER_EVENTMAP_H

#include "Allele.h"
#include "Haplotype.h"
#include "variantcontext/VariantContext.h"

class Haplotype;

class VariantContextComparator;

class EventMap {
private:
    static const int MAX_EVENTS_PER_HAPLOTYPE = 3;
    static const int MAX_INDELS_PER_HAPLOTYPE = 2;
    std::shared_ptr<Haplotype> haplotype;
    std::shared_ptr<uint8_t[]> ref;
    int refLength;
    Locatable* refLoc;
    std::string sourceNameToAdd;
    std::map<int, std::shared_ptr<VariantContext>> variantMap;

public:
    static const Allele* SYMBOLIC_UNASSEMBLED_EVENT_ALLELE;
    EventMap(std::shared_ptr<Haplotype> haplotype, std::shared_ptr<uint8_t[]> ref, int refLength, Locatable* refLoc, std::string sourceNameToAdd, int maxMnpDistance);
    void addVC(std::shared_ptr<VariantContext> vc, bool merge);
    bool empty();
    static std::set<int> buildEventMapsForHaplotypes(std::vector<std::shared_ptr<Haplotype>> & haplotypes, std::shared_ptr<uint8_t[]> ref, int refLength, Locatable* refLoc, bool debug, int maxMnpDistance);
    std::set<int> getStartPositions();
    std::vector<std::shared_ptr<VariantContext>> getVariantContexts();
    static std::set<std::shared_ptr<VariantContext>, VariantContextComparator> getAllVariantContexts(std::vector<std::shared_ptr<Haplotype>> & haplotypes);

protected:
    static const int MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION = 3;
    void processCigarForInitialEvents(int maxMnpDistance);
    std::shared_ptr<VariantContext> makeBlock(std::shared_ptr<VariantContext> vc1, std::shared_ptr<VariantContext> vc2);

};

class VariantContextComparator {
public:
    bool operator() (std::shared_ptr<VariantContext> a1, std::shared_ptr<VariantContext> a2) {
        return a1->getStart() < a2->getStart();
    }
};


#endif //MUTECT2CPP_MASTER_EVENTMAP_H
