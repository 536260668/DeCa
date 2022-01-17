//
// Created by 梦想家xixi on 2021/12/4.
//

#include "EventMap.h"

#include <utility>
#include <vector>
#include "param/ParamUtils.h"
#include "utils/BaseUtils.h"
#include "variantcontext/builder/VariantContextBuilder.h"
#include <deque>

const Allele* EventMap::SYMBOLIC_UNASSEMBLED_EVENT_ALLELE = Allele::create((uint8_t *) "<UNASSEMBLED_EVENT>", 19, false);

EventMap::EventMap(std::shared_ptr<Haplotype> haplotype, uint8_t *ref, int refLength, Locatable *refLoc, std::string sourceNameToAdd,
                   int maxMnpDistance) : haplotype(std::move(haplotype)), ref(ref), refLength(refLength), refLoc(refLoc), sourceNameToAdd(std::move(sourceNameToAdd)){
    processCigarForInitialEvents(maxMnpDistance);
}

void EventMap::processCigarForInitialEvents(int maxMnpDistance) {
    ParamUtils::isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");
    std::shared_ptr<Cigar> cigar = haplotype->getCigar();
    uint8_t * alignment = haplotype->getBases();
    int alignmentLength = haplotype->getLength();

    int refPos = haplotype->getAlignmentStartHapwrtRef();
    if(refPos < 0) {
        return;
    }
    std::vector<std::shared_ptr<VariantContext>> proposedEvents;
    int alignmentPos = 0;
    for(int cigarIndex = 0; cigarIndex < cigar->numCigarElements(); cigarIndex++) {
        CigarElement ce = cigar->getCigarElement(cigarIndex);
        int elementLength = ce.getLength();
        switch (ce.getOperator()) {
            case I:
            {
                if(refPos > 0) {
                    std::vector<Allele*> insertionAlleles;
                    int insertionStart = refLoc->getStart() + refPos - 1;
                    uint8_t refByte = ref[refPos - 1];
                    if(BaseUtils::isRegularBase(refByte)) {
                        insertionAlleles.emplace_back(Allele::create(refByte, true));
                    }
                    if(cigarIndex == 0 || cigarIndex == cigar->numCigarElements() - 1){

                    } else {
                        uint8_t * insertionBases = new uint8_t[1];
                        insertionBases[0] = ref[refPos - 1];
                        int toAddLength;
                        uint8_t * toAdd = Mutect2Utils::copyOfRange(alignment, alignmentLength, alignmentPos, alignmentPos + elementLength, toAddLength);
                        int length = 1+toAddLength;
                        uint8_t * tmp = new uint8_t[length];
                        memcpy(tmp, insertionBases, 1);
                        memcpy(tmp+1, toAdd, toAddLength);
                        delete[] insertionBases;
                        delete[] toAdd;
                        insertionBases = tmp;
                        if(BaseUtils::isAllRegularBases(insertionBases, length)) {
                            insertionAlleles.emplace_back(Allele::create(insertionBases, length, false));
                        }
                    }
                    if(insertionAlleles.size() == 2) {
                        std::string contig = refLoc->getContig();
                        proposedEvents.emplace_back(VariantContextBuilder(sourceNameToAdd, contig, insertionStart, insertionStart, &insertionAlleles).make());
                    }
                }
                alignmentPos += elementLength;
                break;

            }
            case S:
            {
                alignmentPos += elementLength;
                break;
            }
            case D:
            {
                if(refPos > 0) {
                    int deletionBasesLength;
                    uint8_t * deletionBases = Mutect2Utils::copyOfRange(ref, refLength, refPos - 1, refPos + elementLength, deletionBasesLength);
                    std::vector<Allele*> deletionAlleles;
                    int deletionStart = refLoc->getStart() + refPos - 1;
                    uint8_t refByte = ref[refPos - 1];
                    if(BaseUtils::isRegularBase(refByte) && BaseUtils::isAllRegularBases(deletionBases, deletionBasesLength)) {
                        deletionAlleles.emplace_back(Allele::create(deletionBases, deletionBasesLength, true));
                        deletionAlleles.emplace_back(Allele::create(refByte, false));
                        std::string contig = refLoc->getContig();
                        proposedEvents.emplace_back(VariantContextBuilder(sourceNameToAdd, contig, deletionStart, deletionStart + elementLength, &deletionAlleles).make());
                    }
                }
                refPos += elementLength;
                break;
            }
            case M:
            case EQ:
            case X:
            {
                std::deque<int> mismatchOffsets;
                for(int offset = 0; offset < elementLength; offset++) {
                    uint8_t refByte = ref[refPos + offset];
                    uint8_t altByte = alignment[alignmentPos + offset];
                    bool mismatch = refByte != altByte && BaseUtils::isRegularBase(refByte) && BaseUtils::isRegularBase(altByte);
                    if(mismatch) {
                        mismatchOffsets.push_back(offset);
                    }
                }
                while(!mismatchOffsets.empty()) {
                    int start = mismatchOffsets.front();
                    mismatchOffsets.pop_front();
                    int end = start;
                    while (!mismatchOffsets.empty() && mismatchOffsets.front() - end <= maxMnpDistance) {
                        end = mismatchOffsets.front();
                        mismatchOffsets.pop_front();
                    }
                    int length;
                    uint8_t * bases = Mutect2Utils::copyOfRange(ref, refLength, refPos + start, refPos + end + 1, length);
                    Allele* refAllele = Allele::create(bases, length, true);
                    bases = Mutect2Utils::copyOfRange(alignment, alignmentLength, alignmentPos + start, alignmentPos + end + 1, length);
                    Allele* altAllele = Allele::create(bases, length, false);
                    std::string contig = refLoc->getContig();
                    proposedEvents.emplace_back(VariantContextBuilder(sourceNameToAdd, contig, refLoc->getStart() + refPos + start, refLoc->getStart() + refPos + end, new std::vector<Allele*>{refAllele,altAllele}).make());
                }
                refPos += elementLength;
                alignmentPos += elementLength;
                break;
            }
            case N:
            case H:
            case P:
            default:
                throw std::invalid_argument("Unsupported cigar operator created during SW alignment");
        }
    }
    for(std::shared_ptr<VariantContext> proposedEvent : proposedEvents )
        addVC(proposedEvent, true);
}

void EventMap::addVC(std::shared_ptr<VariantContext> vc, bool merge) {
    Mutect2Utils::validateArg(vc.get(), "null is not allowed here");
    if(variantMap.find(vc->getStart()) != variantMap.end()) {
        Mutect2Utils::validateArg(merge, "Will not merge previously bound variant contexts as merge is false");
        std::shared_ptr<VariantContext> prev = variantMap.at(vc->getStart());
        variantMap.insert(std::pair<int, std::shared_ptr<VariantContext>>(vc->getStart(), makeBlock(prev, vc)));
    } else
        variantMap.insert(std::pair<int, std::shared_ptr<VariantContext>>(vc->getStart(), vc));
}

std::shared_ptr<VariantContext> EventMap::makeBlock(std::shared_ptr<VariantContext> vc1, std::shared_ptr<VariantContext> vc2) {
    Mutect2Utils::validateArg(vc1->getStart() == vc2->getStart(), "vc1 and 2 must have the same start");
    Mutect2Utils::validateArg(vc1->isBiallelic(), "vc1 must be biallelic");
    if( ! vc1->isSNP()) {
        Mutect2Utils::validateArg((vc1->isSimpleDeletion() && vc2->isSimpleInsertion()) || (vc1->isSimpleInsertion() && vc2->isSimpleDeletion()), "Can only merge single insertion with deletion (or vice versa)");
    } else {
        Mutect2Utils::validateArg(!vc2->isSNP(), "there's been some terrible bug in the cigar");
    }

    Allele* new_ref, * new_alt;
    VariantContextBuilder b(vc1);
    if(vc1->isSNP()) {
        if(*(vc1->getReference()) == (*(vc2->getReference()))) {
            new_ref = vc1->getReference();
            uint8_t * vc1Bases = vc1->getAlternateAllele(0)->getBases();
            int vc1Length = vc1->getAlternateAllele(0)->getBasesLength();
            uint8_t * vc2Bases = vc2->getAlternateAllele(0)->getBases();
            int vc2Length = vc2->getAlternateAllele(0)->getBasesLength();
            int tmpLength = vc1Length + vc2Length - 1;
            uint8_t * tmp = new uint8_t[tmpLength];
            memcpy(tmp, vc1Bases, vc1Length);
            memcpy(tmp+vc1Length, vc2Bases+1, vc2Length-1);
            new_alt = Allele::create(tmp, tmpLength, false);
        } else {
            new_ref = vc2->getReference();
            new_alt = vc1->getAlternateAllele(0);
            b.setStop(vc2->getEnd());
        }
    } else {
        std::shared_ptr<VariantContext> insertion = vc1->isSimpleInsertion() ? vc1 : vc2;
        std::shared_ptr<VariantContext> deletion  = vc1->isSimpleInsertion() ? vc2 : vc1;
        new_ref = deletion->getReference();
        new_alt = insertion->getAlternateAllele(0);
        b.setStop(deletion->getEnd());
    }
    std::vector<Allele*> alleles{new_alt, new_ref};
    b.setAlleles(&alleles);
    return b.make();
}

bool EventMap::empty() {
    return variantMap.empty();
}

std::set<int> EventMap::buildEventMapsForHaplotypes(std::vector<std::shared_ptr<Haplotype>> & haplotypes, uint8_t *ref, int refLength,
                                                    Locatable *refLoc, bool debug, int maxMnpDistance) {
    Mutect2Utils::validateArg(maxMnpDistance >= 0, "maxMnpDistance may not be negative.");
    std::set<int> startPosKeySet;
    int hapNumber = 0;
    for(const std::shared_ptr<Haplotype>& h : haplotypes) {
        h->setEventMap(new EventMap(h, ref, refLength, refLoc, "HC" + std::to_string(hapNumber), maxMnpDistance));
        for(int i : h->getEventMap()->getStartPositions()) {
            startPosKeySet.insert(i);
        }
    }
    return startPosKeySet;
}

std::set<int> EventMap::getStartPositions() {
    std::set<int> ret;
    for(std::pair<int, std::shared_ptr<VariantContext>> pair_vc : variantMap) {
        ret.insert(pair_vc.first);
    }
    return {ret.begin(), ret.end()};
}

std::set<std::shared_ptr<VariantContext>, VariantContextComparator>
EventMap::getAllVariantContexts(std::vector<std::shared_ptr<Haplotype>> & haplotypes) {
    std::set<std::shared_ptr<VariantContext>, VariantContextComparator> ret;
    for(const std::shared_ptr<Haplotype>& h : haplotypes) {
        for(std::shared_ptr<VariantContext> vc : h->getEventMap()->getVariantContexts()) {
            ret.insert(vc);
        }
    }
    return ret;
}

std::vector<std::shared_ptr<VariantContext>> EventMap::getVariantContexts() {
    std::vector<std::shared_ptr<VariantContext>> ret;
    for(std::pair<int, std::shared_ptr<VariantContext>> pair_vc : variantMap) {
        ret.emplace_back(pair_vc.second);
    }
    return ret;
}
