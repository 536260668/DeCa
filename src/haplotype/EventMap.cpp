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

EventMap::EventMap(Haplotype *haplotype, uint8_t *ref, int refLength, Locatable *refLoc, std::string sourceNameToAdd,
                   int maxMnpDistance) : haplotype(haplotype), ref(ref), refLength(refLength), refLoc(refLoc), sourceNameToAdd(std::move(sourceNameToAdd)){
    processCigarForInitialEvents(maxMnpDistance);
}

void EventMap::processCigarForInitialEvents(int maxMnpDistance) {
    ParamUtils::isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");
    Cigar* cigar = haplotype->getCigar();
    uint8_t * alignment = haplotype->getBases();
    int alignmentLength = haplotype->getLength();

    int refPos = haplotype->getAlignmentStartHapwrtRef();
    if(refPos < 0) {
        return;
    }
    std::vector<VariantContext*> proposedEvents;
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
    for(VariantContext* proposedEvent : proposedEvents )
        addVC(proposedEvent, true);
}

void EventMap::addVC(VariantContext *vc, bool merge) {
    Mutect2Utils::validateArg(vc, "null is not allowed here");
    if(variantMap.find(vc->getStart()) != variantMap.end()) {
        Mutect2Utils::validateArg(merge, "Will not merge previously bound variant contexts as merge is false");
        VariantContext* prev = variantMap.at(vc->getStart());
    }
}

VariantContext *EventMap::makeBlock(VariantContext *vc1, VariantContext *vc2) {
    Mutect2Utils::validateArg(vc1->getStart() == vc2->getStart(), "vc1 and 2 must have the same start");
    Mutect2Utils::validateArg(vc1->isBiallelic(), "vc1 must be biallelic");
    if( ! vc1->isSNP()) {
        Mutect2Utils::validateArg((vc1->isSimpleDeletion() && vc2->isSimpleInsertion()) || (vc1->isSimpleInsertion() && vc2->isSimpleDeletion()), "Can only merge single insertion with deletion (or vice versa)");
    } else {
        Mutect2Utils::validateArg(!vc2->isSNP(), "there's been some terrible bug in the cigar");
    }

    Allele* new_ref, * new_alt;

}
