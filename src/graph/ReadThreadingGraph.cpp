//
// Created by 梦想家xixi on 2021/10/18.
//

#include "ReadThreadingGraph.h"
#include <string>
#include <utility>
#include "Mutect2Utils.h"
#include "smithwaterman/SmithWatermanAlignment.h"
#include "smithwaterman/SWNativeAlignerWrapper.h"
#include "read/AlignmentUtils.h"
#include "SeqVertex.h"

void ReadThreadingGraph::addRead(SAMRecord read) {
    uint8_t* sequence = read.getBases();
    uint8_t* qualities = read.getBaseQualities();

    int lastGood = -1;
    for(int end = 0; end <= read.getLength(); end++) {
        if (end == read.getLength() || !baseIsUsableForAssembly(sequence[end], qualities[end])) {
            int start = lastGood;
            int len = end - start;

            if(start != -1 && len >= kmerSize) {
                std::string name = read.getReadName();
                name += '_' + std::to_string(start) + '_' + std::to_string(end);
                addSequence(name, (std::string &) "SAMFileHeader{VN=1.6, GO=none, SO=coordinate}", sequence, start, end, 1, false);
            }
            lastGood = -1;
        } else if(lastGood == -1) {
            lastGood = end;
        }
    }

}

bool ReadThreadingGraph::baseIsUsableForAssembly(uint8_t base, uint8_t qual) const {
    return base != 'N' && qual >= minBaseQualityToUseInAssembly;
}

void
ReadThreadingGraph::addSequence(std::string seqName, std::string& sampleName, const uint8_t *sequence, int start, int stop,
                                int count, bool isRef) {
    Mutect2Utils::validateArg(!alreadyBuilt, "Graph already built");
    std::map<std::string, std::vector<SequenceForKmers>>::iterator iter;
    iter = pending.find(sampleName);
    if(iter == pending.end()) {
        std::vector<SequenceForKmers> list;
        pending.insert(std::pair<std::string, std::vector<SequenceForKmers>>(sampleName, list));
        iter = pending.find(sampleName);
    }
    iter->second.push_back(SequenceForKmers{std::move(seqName), const_cast<uint8_t *>(sequence), start, stop, count, isRef});
}

std::vector<Kmer> ReadThreadingGraph::determineNonUniqueKmers(SequenceForKmers &sequenceForKmers, const int kmerSize) {
    ArraySet<Kmer> allKmers;
    std::vector<Kmer> nonUniqueKmers;
    const int stopPosition = sequenceForKmers.stop - kmerSize;
    for(int i = 0; i <= stopPosition; i++) {
        Kmer* kmer = new Kmer(sequenceForKmers.sequence, i, kmerSize);
        if(!allKmers.insert(*kmer).second) {
            nonUniqueKmers.push_back(*kmer);
        }
        delete kmer;
    }
    return nonUniqueKmers;
}

std::vector<SequenceForKmers> ReadThreadingGraph::getAllPendingSequences(){
    std::map<std::string, std::vector<SequenceForKmers>>::iterator iter;
    std::vector<SequenceForKmers> res;
    for(iter = pending.begin(); iter != pending.end(); ++iter)
    {
        std::vector<SequenceForKmers>::iterator viter;
        for(viter = iter->second.begin(); viter != iter->second.end(); ++viter)
        {
            res.push_back(*viter);
        }
    }
    return res;
}

ArraySet<Kmer> ReadThreadingGraph::determineKmerSizeAndNonUniques(const int minKmerSize, const int maxKmerSize) {
    std::vector<SequenceForKmers> withNonUniques = getAllPendingSequences();
    ArraySet<Kmer> nonUniqueKmers_m;

    for(int kmerSize_m = minKmerSize; kmerSize_m <= maxKmerSize; kmerSize_m++) {
        nonUniqueKmers_m.clear();
        std::vector<SequenceForKmers>::iterator viter;
        for(viter = withNonUniques.begin(); viter != withNonUniques.end();)
        {
            std::vector<Kmer> nonUniquesFromSeq = determineNonUniqueKmers(*viter, kmerSize_m);
            if(nonUniquesFromSeq.empty()) {
                withNonUniques.erase(viter);
            }
            else {
                for(auto & iter : nonUniquesFromSeq) {
                    nonUniqueKmers_m.insert(iter);
                }
                viter++;
            }
        }

        if(nonUniqueKmers_m.empty())
            break;
    }

    return nonUniqueKmers_m;
}

MultiDeBruijnVertex* ReadThreadingGraph::createVertex(Kmer & kmer) {
    MultiDeBruijnVertex* newVertex = new MultiDeBruijnVertex(kmer.getBases(), kmer.getLength(), false);
    unsigned prevSize = getVertexSet().size();
    addVertex(newVertex);
    if(getVertexSet().size() != prevSize + 1)
        throw std::invalid_argument("Adding vertex to graph didn't increase the graph size");

    if(nonUniqueKmers.find(kmer) == nonUniqueKmers.end() && uniqueKmers.find(kmer) == uniqueKmers.end())
        uniqueKmers.insert(std::pair<Kmer, MultiDeBruijnVertex*>(kmer, newVertex));

    return newVertex;
}

MultiDeBruijnVertex *
ReadThreadingGraph::extendChainByOne(MultiDeBruijnVertex* prevVertex, uint8_t *sequence, const int kmerStart, const int count,
                                     const bool isRef) {
    ArraySet<MultiSampleEdge*> outgoingEdges = outgoingEdgesOf(prevVertex);
    int nextPos = kmerStart + kmerSize - 1;
    std::vector<MultiSampleEdge*>::iterator iter;
    for(iter = outgoingEdges.begin(); iter != outgoingEdges.end(); iter++) {
        MultiDeBruijnVertex* target = getEdgeTarget(*iter);
        if(target->getSuffix() == sequence[nextPos]) {
            (*iter)->incMultiplicity(count);
            return target;
        }
    }

    Kmer *kmer = new Kmer(sequence, kmerStart, kmerSize);
    if(!isRef && *kmer == refSource) {
        MultiDeBruijnVertex* nextVertex = createVertex(*kmer);
        addEdge(prevVertex, nextVertex, new MultiSampleEdge(isRef, count, numPruningSamples));
        delete kmer;
        return nextVertex;
    }
    else {
        if(isRef && getUniqueKmerVertex(*kmer, false))
            throw std::invalid_argument("Found a unique vertex to merge into the reference graph");
        MultiDeBruijnVertex* nextVertex = getUniqueKmerVertex(*kmer, false) ? uniqueKmers.find(*kmer)->second : createVertex(*kmer);
        addEdge(prevVertex, nextVertex, new MultiSampleEdge(isRef, count, numPruningSamples));
        delete kmer;
        return nextVertex;
    }

}

void ReadThreadingGraph::threadSequence(SequenceForKmers & sequenceForKmers) {
    int uniqueStartPos = findStart(sequenceForKmers);
    if(uniqueStartPos == -1)
        return;

    MultiDeBruijnVertex* startingVertex = getOrCreateKmerVertex(sequenceForKmers.sequence, uniqueStartPos);

    if(INCREASE_COUNTS_BACKWARDS) {
        increaseCountsInMatchedKmers(sequenceForKmers, startingVertex, startingVertex->getSequence(), kmerSize - 2);
    }

    if(sequenceForKmers.isRef) {
        if(refSource.getBases() != nullptr)
            throw std::invalid_argument("Found two refSources! prev:");

        refSource = Kmer(sequenceForKmers.sequence, sequenceForKmers.start, kmerSize);
    }

    MultiDeBruijnVertex* vertex = startingVertex;
    for(int i = uniqueStartPos + 1; i <= sequenceForKmers.stop - kmerSize; i++) {
        vertex = extendChainByOne(vertex, sequenceForKmers.sequence, i, sequenceForKmers.count, sequenceForKmers.isRef);
    }
}

int ReadThreadingGraph::findStart(SequenceForKmers seqForKmers) {
    if(seqForKmers.isRef)
        return 0;

    for(int i = seqForKmers.start; i < seqForKmers.stop - kmerSize; i++) {
        Kmer kmer1(seqForKmers.sequence, i, kmerSize);
        if((startThreadingOnlyAtExistingVertex ? uniqueKmers.find(kmer1) != uniqueKmers.end() : nonUniqueKmers.find(kmer1) == nonUniqueKmers.end()))
            return i;
    }

    return -1;
}

bool ReadThreadingGraph::getUniqueKmerVertex(Kmer & kmer, const bool allowRefSource) {
    if(!allowRefSource && kmer == refSource)
        return false;
    else {
        bool res = uniqueKmers.find(kmer) != uniqueKmers.end();
        return res;
    }


}

MultiDeBruijnVertex* ReadThreadingGraph::getOrCreateKmerVertex(uint8_t *sequence, const int start) {
    Kmer kmer(sequence, start, kmerSize);
    return getUniqueKmerVertex(kmer, true) ? uniqueKmers.find(kmer)->second : createVertex(kmer);
}

void ReadThreadingGraph::increaseCountsInMatchedKmers(SequenceForKmers & seqForKmers, MultiDeBruijnVertex* vertex,
                                                      uint8_t* originalKmer, int offset) {
    if(offset == -1)
        return;

    ArraySet<MultiSampleEdge*> incomingEdges = incomingEdgesOf(vertex);
    for(std::vector<MultiSampleEdge*>::iterator iter = incomingEdges.begin(); iter != incomingEdges.end(); iter++) {
        MultiDeBruijnVertex* prev = getEdgeSource(*iter);
        uint8_t suffix = prev->getSuffix();
        uint8_t seqBase = originalKmer[offset];
        if(suffix == seqBase && (increaseCountsThroughBranches || inDegreeOf(vertex) == 1)) {
            (*iter)->incMultiplicity(seqForKmers.count);
            increaseCountsInMatchedKmers(seqForKmers, prev, originalKmer, offset-1);
        }
    }
}

void ReadThreadingGraph::buildGraphIfNecessary() {
    if(alreadyBuilt)
        return;

    nonUniqueKmers = determineKmerSizeAndNonUniques(kmerSize, kmerSize);


    for(std::map<std::string, std::vector<SequenceForKmers>>::iterator miter = pending.begin(); miter != pending.end(); miter++){
        for(std::vector<SequenceForKmers>::iterator viter = miter->second.begin(); viter != miter->second.end(); viter++) {
            threadSequence(*viter);
        }
        std::map<MultiSampleEdge*, IntrusiveEdge<MultiDeBruijnVertex>>::iterator eiter;
        for(eiter = edgeMap.begin(); eiter != edgeMap.end(); eiter++) {
            (*eiter->first).flushSingleSampleMultiplicity();
        }
    }


    pending.clear();
    alreadyBuilt = true;

    for(std::map<Kmer, MultiDeBruijnVertex*>::iterator kiter = uniqueKmers.begin(); kiter != uniqueKmers.end(); kiter++) {
        kiter->second->setAdditionalInfo(kiter->second->getAdditionalInfo() + '+');
    }

}

void ReadThreadingGraph::setThreadingStartOnlyAtExistingVertex(bool value) {
    startThreadingOnlyAtExistingVertex = value;
}

void ReadThreadingGraph::setPending() {
    uint8_t* byte997137 = new uint8_t[99]{67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84};
    SequenceForKmers sequenceForKmers7135 = {.name = "HWI-ST729_110151799:2:41:9503:140222_0_99", .sequence = byte997137, .start = 0, .stop = 99, .count = 1, .isRef = false};
    std::vector<SequenceForKmers> v3;
    v3.emplace_back(sequenceForKmers7135);
    std::string key1491293244 = "H_LS-E2-A15C-01A-31D-A12B-09";

    uint8_t* byte1007126 = new uint8_t[100]{67, 84, 65, 65, 67, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67};
    SequenceForKmers sequenceForKmers7123 = {.name = "HWI-ST729_110151799:3:66:3953:197183_0_100", .sequence = byte1007126, .start = 0, .stop = 100, .count = 1, .isRef = false};
    uint8_t* byte907131 = new uint8_t[90]{65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 67, 67, 65};
    SequenceForKmers sequenceForKmers7124 = {.name = "HWI-ST729_110151799:3:64:7208:158288_0_90", .sequence = byte907131, .start = 0, .stop = 90, .count = 1, .isRef = false};
    std::vector<SequenceForKmers> v2;
    v2.emplace_back(sequenceForKmers7123);
    v2.emplace_back(sequenceForKmers7124);
    std::string key1617703971 = "H_LS-E2-A15C-10A-01D-A12B-09";

    uint8_t* byte2867095 = new uint8_t[286]{67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 65, 67, 67, 67, 84, 65, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67};
    SequenceForKmers sequenceForKmers7094 = {.name = "ref", .sequence = byte2867095, .start = 0, .stop = 286, .count = 1, .isRef = true};
    std::vector<SequenceForKmers> v1;
    v1.emplace_back(sequenceForKmers7094);
    std::string key769618002 = "XXX_UNNAMED_XXX";
    pending.insert(std::pair<std::string, std::vector<SequenceForKmers>>(key1617703971, v2));
    pending.insert(std::pair<std::string, std::vector<SequenceForKmers>>(key769618002, v3));
    pending.insert(std::pair<std::string, std::vector<SequenceForKmers>>(key1491293244, v1));
}

bool ReadThreadingGraph::removeVertex(MultiDeBruijnVertex *V) {
    uint8_t * sequence = new uint8_t[V->getLength()];
    memcpy(sequence, V->getSequence(), V->getLength());
    bool result = DirectedSpecifics::removeVertex(V);

    if(result) {
        Kmer kmer(sequence, 0 ,kmerSize);
        uniqueKmers.erase(kmer);
    }
    delete[] sequence;

    return result;
}

void ReadThreadingGraph::removeSingletonOrphanVertices() {
    std::vector<MultiDeBruijnVertex*> toRemove;
    ArraySet<MultiDeBruijnVertex*> allvertex = getVertexSet();
    typename std::vector<MultiDeBruijnVertex*>::iterator viter;
    for(viter = allvertex.begin(); viter != allvertex.end(); viter++) {
        if(inDegreeOf(*viter) == 0 && outDegreeOf(*viter) == 0) {
            toRemove.emplace_back(*viter);
        }
    }
    removeAllVertices(toRemove);
}

void ReadThreadingGraph::recoverDanglingTails(int pruneFactor, int minDanglingBranchLength, bool recoverAll) {
    Mutect2Utils::validateArg(pruneFactor >= 0, "pruneFactor must be non-negative");
    Mutect2Utils::validateArg(minDanglingBranchLength >= 0, "minDanglingBranchLength must be non-negative");

    if(!alreadyBuilt) {
        throw std::invalid_argument("recoverDanglingTails requires the graph be already built");
    }

    int attempted = 0;
    int nRecovered = 0;
    for ( MultiDeBruijnVertex* v :  getVertexSet()) {
        if(outDegreeOf(v) == 0 && !isRefSink(v)) {
            attempted++;
            nRecovered += recoverDanglingTail(v, pruneFactor, minDanglingBranchLength, recoverAll);
        }
    }
}

int ReadThreadingGraph::recoverDanglingTail(MultiDeBruijnVertex *vertex, int pruneFactor, int minDanglingBranchLength,
                                            bool recoverAll) {
    if(outDegreeOf(vertex) != 0 ) {
        throw std::invalid_argument("Attempting to recover a dangling tail but it has out-degree > 0");
    }

    DanglingChainMergeHelper* danglingTailMergeResult = generateCigarAgainstDownwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength, recoverAll);

    if(danglingTailMergeResult == nullptr || !cigarIsOkayToMerge(danglingTailMergeResult->cigar, false, true)) {
        return 0;
    }
    return mergeDanglingTail(danglingTailMergeResult);
}

DanglingChainMergeHelper*
ReadThreadingGraph::generateCigarAgainstDownwardsReferencePath(MultiDeBruijnVertex *vertex, int pruneFactor,
                                                               int minDanglingBranchLength, bool recoverAll) {
    int minTailPathLength = std::max(1, minDanglingBranchLength);
    std::deque<MultiDeBruijnVertex*> altPath = findPathUpwardsToLowestCommonAncestor(vertex, pruneFactor, !recoverAll);
    if(altPath.empty() || isRefSource(altPath.front()) || altPath.size() < minTailPathLength + 1) {
        return nullptr;
    }
    MultiDeBruijnVertex* toBlacklistedEdge =  altPath.size() > 1 ? altPath.at(1) : nullptr;
    std::vector<MultiDeBruijnVertex*> refPath = getReferencePath(altPath.at(0), downwards, getHeaviestIncomingEdge(toBlacklistedEdge));
    std::vector<MultiDeBruijnVertex*> newaltPath;
    for(MultiDeBruijnVertex* vertex1 : altPath) {
        newaltPath.emplace_back(vertex1);
    }
    int refLength;
    uint8_t * refBases = getBasesForPath(refPath, refLength, false);
    int altLength;
    uint8_t * altBases = getBasesForPath(newaltPath, altLength, false);
    SWNativeAlignerWrapper wrapper = SWNativeAlignerWrapper();
    SmithWatermanAlignment* alignment = wrapper.align(refBases, refLength, altBases, altLength, &SmithWatermanAligner::STANDARD_NGS, LEADING_INDEL);
    return new DanglingChainMergeHelper(newaltPath, refPath, altBases, altLength, refBases, refLength, AlignmentUtils::removeTrailingDeletions(alignment->getCigar()));
}



std::deque<MultiDeBruijnVertex *>
ReadThreadingGraph::findPathUpwardsToLowestCommonAncestor(MultiDeBruijnVertex *vertex, int pruneFactor,
                                                          bool giveUpAtBranch) {
    std::deque<MultiDeBruijnVertex*> ret;
    MultiDeBruijnVertex* v = vertex;
    if(giveUpAtBranch) {
        while(!(inDegreeOf(v) != 1 || outDegreeOf(v) >= 2)) {
            MultiSampleEdge* edge = incomingEdgeOf(v);
            if(edge->getPruningMultiplicity() < pruneFactor) {
                ret.clear();
            }
            else {
                ret.push_front(v);
            }
            v = getEdgeSource(edge);
        }
        ret.push_front(v);

        return outDegreeOf(v) > 1 ? ret : std::deque<MultiDeBruijnVertex*>();
    } else {
        while(!(hasIncidentRefEdge(v) || inDegreeOf(v) == 0)) {
            MultiSampleEdge* edge = getHeaviestIncomingEdge(v);
            if(edge->getPruningMultiplicity() < pruneFactor) {
                ret.clear();
            }
            else {
                ret.push_front(v);
            }
            v = getEdgeSource(edge);
        }
        ret.push_front(v);

        return outDegreeOf(v) > 1 && hasIncidentRefEdge(v) ? ret : std::deque<MultiDeBruijnVertex*>();
    }
}

bool ReadThreadingGraph::hasIncidentRefEdge(MultiDeBruijnVertex *v) {
    for(MultiSampleEdge* edge : incomingEdgesOf(v)) {
        if(edge->getIsRef()) {
            return true;
        }
    }
    return false;
}

MultiSampleEdge *ReadThreadingGraph::getHeaviestIncomingEdge(MultiDeBruijnVertex *v) {
    ArraySet<MultiSampleEdge*> incomingEdges = incomingEdgesOf(v);
    MultiSampleEdge* ret;
    ret = *incomingEdges.begin();
    for(std::vector<MultiSampleEdge*>::iterator iter = incomingEdges.begin(); iter != incomingEdges.end(); iter++) {
        if(ret->getPruningMultiplicity() < (*iter)->getPruningMultiplicity())
            ret = *iter;
    }
    return ret;
}

std::vector<MultiDeBruijnVertex *>
ReadThreadingGraph::getReferencePath(MultiDeBruijnVertex *start, TraversalDirection direction,
                                     MultiSampleEdge *blacklistedEdge) {
    std::vector<MultiDeBruijnVertex *> path;

    MultiDeBruijnVertex* v = start;

    while (v != nullptr) {
        path.emplace_back(v);
        v = (direction == downwards ? getNextReferenceVertex(v, true, blacklistedEdge) : getPrevReferenceVertex(v));
    }
    return path;
}

uint8_t *ReadThreadingGraph::getBasesForPath(std::vector<MultiDeBruijnVertex *> path, int &length, bool expandSource) {
    int tmpLength = 300;
    int start = 0;
    uint8_t * tmp = new uint8_t[300];
    for(MultiDeBruijnVertex* v : path) {
        if(expandSource && isSource(v)) {
            uint8_t * seq = v->getSequence();
            int seqLength = v->getLength();
            while(start + seqLength > tmpLength) {
                tmpLength *= 2;
                uint8_t * newtmp = new uint8_t[tmpLength];
                memcpy(newtmp, tmp, tmpLength);
                delete[] tmp;
                tmp = newtmp;
            }
            for(int i = 0; i < seqLength; i++) {
                tmp[start + i] = seq[seqLength-1-i];
            }
            start += seqLength;
        } else {
            while(start + 1 > tmpLength) {
                tmpLength *= 2;
                uint8_t * newtmp = new uint8_t[tmpLength];
                memcpy(newtmp, tmp, tmpLength);
                delete[] tmp;
                tmp = newtmp;
            }
            tmp[start] = v->getSuffix();
            start++;
        }
    }
    uint8_t* ret = new uint8_t[start+1];
    memcpy(ret, tmp, start);
    ret[start] = '\0';
    delete[] tmp;
    length = start;
    return ret;
}

bool ReadThreadingGraph::cigarIsOkayToMerge(Cigar *cigar, bool requireFirstElementM, bool requireLastElementM) {
    std::vector<CigarElement> elements = cigar->getCigarElements();
    int numElements = elements.size();
    if ( numElements == 0 || numElements > MAX_CIGAR_COMPLEXITY ) {
        return false;
    }

    if ( requireFirstElementM && elements.at(0).getOperator() != M ) {
        return false;
    }

    if ( requireLastElementM && elements.at(numElements - 1).getOperator() != M ) {
        return false;
    }

    return true;
}

int ReadThreadingGraph::mergeDanglingTail(DanglingChainMergeHelper *danglingTailMergeResult) {
    std::vector<CigarElement> elements = danglingTailMergeResult->cigar->getCigarElements();
    CigarElement lastElement = elements.at(elements.size() - 1);
    Mutect2Utils::validateArg(lastElement.getOperator() == M, "The last Cigar element must be an M");

    int lastRefIndex = danglingTailMergeResult->cigar->getReferenceLength() - 1;
    int matchingSuffix = std::min(longestSuffixMatch(danglingTailMergeResult->referencePathString, danglingTailMergeResult->referencePathStringLength, danglingTailMergeResult->danglingPathString, danglingTailMergeResult->danglingPathStringLength, lastRefIndex),
            lastElement.getLength());

    if(matchingSuffix == 0) {
        return 0;
    }

    int altIndexToMerge = std::max(danglingTailMergeResult->cigar->getReadLength() - matchingSuffix - 1, 0);

    bool firstElementIsDeletion = elements.at(0).getOperator() == D;
    bool mustHandleLeadingDeletionCase = firstElementIsDeletion && (elements.at(0).getLength() + matchingSuffix == lastRefIndex + 1);
    int refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (mustHandleLeadingDeletionCase ? 1 : 0);

    if ( refIndexToMerge == 0 ) {
        return 0;
    }
    addEdge(danglingTailMergeResult->danglingPath.at(altIndexToMerge), danglingTailMergeResult->referencePath.at(refIndexToMerge), new MultiSampleEdge(
            false, 1, numPruningSamples));
    return 1;
}

int ReadThreadingGraph::longestSuffixMatch(uint8_t *seq, int seqLength, uint8_t *kmer, int kmerLength, int seqStart) {
    for (int len = 1; len <= kmerLength; len++) {
        int seqI = seqStart - len + 1;
        int kmerI = kmerLength - len;
        if ( seqI < 0 || seq[seqI] != kmer[kmerI] ) {
            return len - 1;
        }
    }
    return kmerLength;
}

void ReadThreadingGraph::recoverDanglingHeads(int pruneFactor, int minDanglingBranchLength, bool recoverAll) {
    Mutect2Utils::validateArg(pruneFactor >= 0, "pruneFactor must be non-negative");
    Mutect2Utils::validateArg(minDanglingBranchLength >= 0, "minDanglingBranchLength must be non-negative");

    if(!alreadyBuilt) {
        throw std::invalid_argument("recoverDanglingTails requires the graph be already built");
    }

    std::vector<MultiDeBruijnVertex*> danglingHeads;
    for ( MultiDeBruijnVertex* v :  getVertexSet()) {
        if( DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::inDegreeOf(v) == 0 && !isRefSource(v))
            danglingHeads.emplace_back(v);
    }

    int attempted = 0;
    int nRecovered = 0;
    for ( MultiDeBruijnVertex* v :  danglingHeads) {
        if(outDegreeOf(v) == 0 && !isRefSink(v)) {
            attempted++;
            nRecovered += recoverDanglingHead(v, pruneFactor, minDanglingBranchLength, recoverAll);
        }
    }
}

int ReadThreadingGraph::recoverDanglingHead(MultiDeBruijnVertex *vertex, int pruneFactor, int minDanglingBranchLength,
                                            bool recoverAll) {
    if(inDegreeOf(vertex) != 0 ) {
        throw std::invalid_argument("Attempting to recover a dangling head but it has in-degree > 0");
    }

    DanglingChainMergeHelper* danglingTailMergeResult = generateCigarAgainstUpwardsReferencePath(vertex, pruneFactor, minDanglingBranchLength, recoverAll);

    if(danglingTailMergeResult == nullptr || !cigarIsOkayToMerge(danglingTailMergeResult->cigar, true, false)) {
        return 0;
    }
    return mergeDanglingHead(danglingTailMergeResult);
}

DanglingChainMergeHelper *
ReadThreadingGraph::generateCigarAgainstUpwardsReferencePath(MultiDeBruijnVertex *vertex, int pruneFactor,
                                                             int minDanglingBranchLength, bool recoverAll) {
    std::deque<MultiDeBruijnVertex*> altPath = findPathDownwardsToHighestCommonDescendantOfReference(vertex, pruneFactor, !recoverAll);
    if(altPath.empty() || isRefSink(altPath.front()) || altPath.size() < minDanglingBranchLength + 1) {
        return nullptr;
    }
    std::vector<MultiDeBruijnVertex*> refPath = getReferencePath(altPath.at(0), upwards, nullptr);
    std::vector<MultiDeBruijnVertex*> newaltPath;
    for(MultiDeBruijnVertex* vertex1 : altPath) {
        newaltPath.emplace_back(vertex1);
    }
    int refLength;
    uint8_t * refBases = getBasesForPath(refPath, refLength, false);
    int altLength;
    uint8_t * altBases = getBasesForPath(newaltPath, altLength, false);
    SWNativeAlignerWrapper wrapper = SWNativeAlignerWrapper();
    SmithWatermanAlignment* alignment = wrapper.align(refBases, refLength, altBases, altLength, &SmithWatermanAligner::STANDARD_NGS, LEADING_INDEL);
    return new DanglingChainMergeHelper(newaltPath, refPath, altBases, altLength, refBases, refLength, AlignmentUtils::removeTrailingDeletions(alignment->getCigar()));
}

std::deque<MultiDeBruijnVertex *>
ReadThreadingGraph::findPathDownwardsToHighestCommonDescendantOfReference(MultiDeBruijnVertex *vertex, int pruneFactor,
                                                                          bool giveUpAtBranch) {
    std::deque<MultiDeBruijnVertex*> ret;
    MultiDeBruijnVertex* v = vertex;
    if(giveUpAtBranch) {
        while(isReferenceNode(v) || outDegreeOf(v) != 1) {
            MultiSampleEdge* edge = outgoingEdgeOf(v);
            if(edge->getPruningMultiplicity() < pruneFactor) {
                ret.clear();
            }
            else {
                ret.push_front(v);
            }
            v = getEdgeTarget(edge);
        }
        ret.push_front(v);

        return isReferenceNode(v)? ret : std::deque<MultiDeBruijnVertex*>();
    } else {
        while(isReferenceNode(v) || outDegreeOf(v) == 0) {
            MultiSampleEdge* edge = getHeaviestOutgoingEdge(v);
            if(edge->getPruningMultiplicity() < pruneFactor) {
                ret.clear();
            }
            else {
                ret.push_front(v);
            }
            v = getEdgeTarget(edge);
        }
        ret.push_front(v);

        return isReferenceNode(v) ? ret : std::deque<MultiDeBruijnVertex*>();
    }
}

MultiSampleEdge *ReadThreadingGraph::getHeaviestOutgoingEdge(MultiDeBruijnVertex *v) {
    ArraySet<MultiSampleEdge*> outgoing = outgoingEdgesOf(v);
    MultiSampleEdge* ret;
    ret = *outgoing.begin();
    for(std::vector<MultiSampleEdge*>::iterator iter = outgoing.begin(); iter != outgoing.end(); iter++) {
        if(ret->getPruningMultiplicity() < (*iter)->getPruningMultiplicity())
            ret = *iter;
    }
    return ret;
}

int ReadThreadingGraph::mergeDanglingHead(DanglingChainMergeHelper *danglingHeadMergeResult) {
    std::vector<CigarElement> elements = danglingHeadMergeResult->cigar->getCigarElements();
    CigarElement firstElement = elements.at(0);
    Mutect2Utils::validateArg(firstElement.getOperator() == M, "The first Cigar element must be an M");
    int indexesToMerge = bestPrefixMatch(danglingHeadMergeResult->referencePathString, danglingHeadMergeResult->danglingPathString, firstElement.getLength());
    if(indexesToMerge <= 0) {
        return 0;
    }

    if(indexesToMerge >= danglingHeadMergeResult->referencePath.size() - 1) {
        return 0;
    }

    if(indexesToMerge >= danglingHeadMergeResult->danglingPath.size() && ! extendDanglingPathAgainstReference(danglingHeadMergeResult, indexesToMerge - danglingHeadMergeResult->danglingPath.size() + 2))
        return 0;

    addEdge(danglingHeadMergeResult->referencePath.at(indexesToMerge+1), danglingHeadMergeResult->danglingPath.at(indexesToMerge), new MultiSampleEdge(
            false, 1, numPruningSamples));

    return 1;
}

int ReadThreadingGraph::bestPrefixMatch(const uint8_t *path1, const uint8_t *path2, int maxIndex) {
    int maxMismatches = getMaxMismatches(maxIndex);
    int mismatches = 0;
    int index = 0;
    int lastGoodIndex = -1;
    while ( index < maxIndex ) {
        if ( path1[index] != path2[index] ) {
            if ( ++mismatches > maxMismatches ) {
                return -1;
            }
            lastGoodIndex = index;
        }
        index++;
    }
    // if we got here then we hit the max index
    return lastGoodIndex;
}

int ReadThreadingGraph::getMaxMismatches(int lengthOfDanglingBranch) const {
    return maxMismatchesInDanglingHead > 0 ? maxMismatchesInDanglingHead : std::max(1, (lengthOfDanglingBranch / kmerSize));
}

bool ReadThreadingGraph::extendDanglingPathAgainstReference(DanglingChainMergeHelper *danglingHeadMergeResult,
                                                            int numNodesToExtend) {
    int indexOfLastDanglingNode = danglingHeadMergeResult->danglingPath.size() - 1;
    int indexOfRefNodeToUse = indexOfLastDanglingNode + numNodesToExtend;
    if ( indexOfRefNodeToUse >= danglingHeadMergeResult->referencePath.size() ) {
        return false;
    }

    MultiDeBruijnVertex* danglingSource = danglingHeadMergeResult->danglingPath.at(indexOfLastDanglingNode);
    danglingHeadMergeResult->danglingPath.erase(danglingHeadMergeResult->danglingPath.begin() + indexOfLastDanglingNode);
    uint8_t * danglingSourceSeq = danglingSource->getSequence();
    int danglingSourceSeqLength = danglingSource->getLength();
    uint8_t * refSourceSequence = danglingHeadMergeResult->referencePath.at(indexOfRefNodeToUse)->getSequence();
    uint8_t * sequenceToExtend = new uint8_t[numNodesToExtend + danglingSourceSeqLength];
    for ( int i = 0; i < numNodesToExtend; i++ ) {
        sequenceToExtend[i] = refSourceSequence[i];
    }
    memcpy(sequenceToExtend+numNodesToExtend, danglingSourceSeq, danglingSourceSeqLength);

    MultiSampleEdge* sourceEdge = getHeaviestOutgoingEdge(danglingSource);
    MultiDeBruijnVertex* prevV = getEdgeTarget(sourceEdge);
    MultiSampleEdge* ret = removeEdge(danglingSource, prevV);
    delete ret;
    for( int i = numNodesToExtend; i > 0; i-- ) {
        uint8_t * tmp = new uint8_t[kmerSize];
        memcpy(tmp, sequenceToExtend + i, kmerSize);
        MultiDeBruijnVertex* newV = new MultiDeBruijnVertex(tmp, kmerSize);
        addVertex(newV);
        MultiSampleEdge* newE = addEdge(newV, prevV);
        newE->setMultiplicity(sourceEdge->getMultiplicity());
        danglingHeadMergeResult->danglingPath.emplace_back(newV);
        prevV = newV;
    }
    return true;
}

MultiSampleEdge *ReadThreadingGraph::createEdge(MultiDeBruijnVertex *, MultiDeBruijnVertex *) {
    return new MultiSampleEdge(false, 1, numPruningSamples);
}

SeqGraph *ReadThreadingGraph::toSequenceGraph() {
    buildGraphIfNecessary();
    SeqGraph* seqGraph = new SeqGraph(kmerSize);
    std::map<MultiDeBruijnVertex*, SeqVertex*> vertexMap;
    for(MultiDeBruijnVertex* dv : DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::getVertexSet()) {
        SeqVertex* sv = new SeqVertex(dv->getAdditionalSequence(
                DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::isSource(dv)), dv->getAdditionalLength(DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::isSource(dv)));
        sv->setAdditionalInfo(dv->getAdditionalInfo());
        vertexMap.insert(std::pair<MultiDeBruijnVertex*, SeqVertex*>(dv, sv));
        seqGraph->addVertex(sv);
    }

    std::map<MultiSampleEdge*, IntrusiveEdge<MultiDeBruijnVertex>>::iterator eiter;
    for(eiter = edgeMap.begin(); eiter != edgeMap.end(); eiter++) {
        SeqVertex* seqInV = vertexMap.at(getEdgeSource(eiter->first));
        SeqVertex* seqOutV = vertexMap.at(getEdgeTarget(eiter->first));
        seqGraph->addEdge(seqInV, seqOutV, new BaseEdge(eiter->first->getIsRef(), eiter->first->getMultiplicity()));
    }
    return seqGraph;
}







