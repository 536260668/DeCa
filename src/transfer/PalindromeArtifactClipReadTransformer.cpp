//
// Created by 梦想家xixi on 2022/2/22.
//

#include "PalindromeArtifactClipReadTransformer.h"

#include <utility>
#include "read/ReadUtils.h"
#include "utils/BaseUtils.h"
#include "clipping/ReadClipper.h"

const double PalindromeArtifactClipReadTransformer::MIN_FRACTION_OF_MATCHING_BASES = 0.9;

PalindromeArtifactClipReadTransformer::PalindromeArtifactClipReadTransformer(
        std::shared_ptr<ReferenceCache> referenceDataSource, SAMFileHeader *header, int minPalindromeSize) : referenceDataSource(std::move(referenceDataSource)), header(header), minPalindromeSize(minPalindromeSize){

}

std::shared_ptr<SAMRecord> PalindromeArtifactClipReadTransformer::apply(const std::shared_ptr<SAMRecord> &read) {
    int adaptorBoundary = read->getAdaptorBoundary();
    if(!read->isProperlyPaired() || adaptorBoundary == ReadUtils::CANNOT_COMPUTE_ADAPTOR_BOUNDARY) {
        return read;
    }
    std::shared_ptr<Cigar> cigar = read->getCigar();
    CigarOperator firstOperator = cigar->getFirstCigarElement().getOperator();
    CigarOperator lastOperator = cigar->getLastCigarElement().getOperator();
    bool readIsUpstreamOfMate = read->getFragmentLength() > 0;

    if((readIsUpstreamOfMate && firstOperator != S && firstOperator != I) ||
            (!readIsUpstreamOfMate && lastOperator != S && lastOperator != I)) {
        return read;
    }

    int potentialArtifactBaseCount = readIsUpstreamOfMate ? cigar->getFirstCigarElement().getLength() :
                                     cigar->getLastCigarElement().getLength();

    int numBasesToCompare = std::min(potentialArtifactBaseCount + minPalindromeSize, read->getLength());

    int contig = header->getSequenceIndex(read->getContig());
    int refStart = readIsUpstreamOfMate ? adaptorBoundary - numBasesToCompare : adaptorBoundary + 1;
    int refEnd = readIsUpstreamOfMate ? adaptorBoundary - 1 : adaptorBoundary + numBasesToCompare;

    if(refStart < 1 || refEnd > header->getSequenceDictionary().getSequences()[contig].getSequenceLength()) {
        return read;
    }

    if((readIsUpstreamOfMate && refStart < read->getStart()) || (!readIsUpstreamOfMate && read->getEnd() < refEnd)) {
        return read;
    }


    int numMatch = 0;

    int readIndex = readIsUpstreamOfMate ? numBasesToCompare - 1 : read->getLength() - 1;
    int length = 0;
    std::shared_ptr<uint8_t[]> refBases = referenceDataSource->getSubsequenceAt(contig, refStart, refEnd, length);

    for(int i = 0; i < length; i++) {
        if(BaseUtils::getComplement(refBases[i]) == read->getBase(readIndex)) {
            numMatch++;
        }
        readIndex--;
    }

    if(numMatch / static_cast<double>(numBasesToCompare) >= MIN_FRACTION_OF_MATCHING_BASES) {
        ReadClipper readClipper(read);
        ClippingOp clippingOp = readIsUpstreamOfMate ? ClippingOp(0, potentialArtifactBaseCount - 1) :
                ClippingOp(read->getLength() - potentialArtifactBaseCount, read->getLength());
        readClipper.addOp(clippingOp);
        return readClipper.clipRead(HARDCLIP_BASES);
    } else {
        return read;
    }
}
