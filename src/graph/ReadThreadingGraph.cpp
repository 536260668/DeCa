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

void ReadThreadingGraph::addRead(std::shared_ptr<SAMRecord> & read) {
    std::shared_ptr<uint8_t[]> sequence_ = read->getBasesNoCopy();
    std::shared_ptr<uint8_t[]> qualities_ = read->getBaseQualitiesNoCopy();
    uint8_t * sequence = sequence_.get();
    uint8_t * qualities = qualities_.get();

    int lastGood = -1;
    for(int end = 0; end <= read->getLength(); end++) {
        if (end == read->getLength() || !baseIsUsableForAssembly(sequence[end], qualities[end])) {
            int start = lastGood;
            int len = end - start;

            if(start != -1 && len >= kmerSize) {
                std::string name = read->getName();
                name += '_' + std::to_string(start) + '_' + std::to_string(end);
                std::string sampleName = read->getGroup() == 0 ? "normal" : "tumor";
                addSequence(name, sampleName, sequence_, start, end, 1, false);
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
ReadThreadingGraph::addSequence(std::string seqName, std::string& sampleName, const std::shared_ptr<uint8_t[]>& sequence, int start, int stop,
                                int count, bool isRef) {
    Mutect2Utils::validateArg(!alreadyBuilt, "Graph already built");
    std::map<std::string, std::vector<SequenceForKmers>>::iterator iter;
    iter = pending.find(sampleName);
    if(iter == pending.end()) {
        std::vector<SequenceForKmers> list;
        pending.insert(std::pair<std::string, std::vector<SequenceForKmers>>(sampleName, list));
        iter = pending.find(sampleName);
    }
    iter->second.push_back(SequenceForKmers{std::move(seqName), sequence, start, stop, count, isRef});
}

std::vector<Kmer> ReadThreadingGraph::determineNonUniqueKmers(SequenceForKmers &sequenceForKmers, const int kmerSize) {
    std::unordered_set<Kmer, hash_kmer, equal_kmer> allKmers;
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

std::list<SequenceForKmers> ReadThreadingGraph::getAllPendingSequences(){
    std::map<std::string, std::vector<SequenceForKmers>>::iterator iter;
    std::list<SequenceForKmers> res;
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

std::unordered_set<Kmer, hash_kmer, equal_kmer> ReadThreadingGraph::determineKmerSizeAndNonUniques(const int minKmerSize, const int maxKmerSize) {
    std::list<SequenceForKmers> withNonUniques = getAllPendingSequences();
    std::unordered_set<Kmer, hash_kmer, equal_kmer> nonUniqueKmers_m;

    for(int kmerSize_m = minKmerSize; kmerSize_m <= maxKmerSize; kmerSize_m++) {
        nonUniqueKmers_m.clear();
        std::list<SequenceForKmers>::iterator viter;
        for(viter = withNonUniques.begin(); viter != withNonUniques.end();)
        {
            std::vector<Kmer> nonUniquesFromSeq = determineNonUniqueKmers(*viter, kmerSize_m);
            if(nonUniquesFromSeq.empty()) {
                withNonUniques.erase(viter++);
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

std::shared_ptr<MultiDeBruijnVertex> ReadThreadingGraph::createVertex(Kmer & kmer) {
    std::shared_ptr<MultiDeBruijnVertex> newVertex(new MultiDeBruijnVertex(kmer.getBases(), kmer.getLength(), false));
    unsigned prevSize = getVertexSet().size();
    addVertex(newVertex);
    if(getVertexSet().size() != prevSize + 1)
        throw std::invalid_argument("Adding vertex to graph didn't increase the graph size");

    if(nonUniqueKmers.find(kmer) == nonUniqueKmers.end() && uniqueKmers.find(kmer) == uniqueKmers.end())
        uniqueKmers.insert(std::pair<Kmer, std::shared_ptr<MultiDeBruijnVertex>>(kmer, newVertex));

    return newVertex;
}

std::shared_ptr<MultiDeBruijnVertex>
ReadThreadingGraph::extendChainByOne(const std::shared_ptr<MultiDeBruijnVertex>& prevVertex, std::shared_ptr<uint8_t[]>sequence, const int kmerStart, const int count,
                                     const bool isRef) {
    std::unordered_set<std::shared_ptr<MultiSampleEdge>> outgoingEdges = outgoingEdgesOf(prevVertex);
    int nextPos = kmerStart + kmerSize - 1;
    std::unordered_set<std::shared_ptr<MultiSampleEdge>>::iterator iter;
    uint8_t * sequence_ = sequence.get();
    for(iter = outgoingEdges.begin(); iter != outgoingEdges.end(); iter++) {
        std::shared_ptr<MultiDeBruijnVertex> target = getEdgeTarget(*iter);
        if(target->getSuffix() == sequence_[nextPos]) {
            (*iter)->incMultiplicity(count);
            return target;
        }
    }

    Kmer *kmer = new Kmer(sequence, kmerStart, kmerSize);
    if(!isRef && *kmer == refSource) {
        std::shared_ptr<MultiDeBruijnVertex> nextVertex = createVertex(*kmer);
        addEdge(prevVertex, nextVertex, std::shared_ptr<MultiSampleEdge>(new MultiSampleEdge(isRef, count, numPruningSamples)));
        delete kmer;
        return nextVertex;
    }
    else {
        if(isRef && getUniqueKmerVertex(*kmer, false))
            throw std::invalid_argument("Found a unique vertex to merge into the reference graph");
        std::shared_ptr<MultiDeBruijnVertex> nextVertex = getUniqueKmerVertex(*kmer, false) ? uniqueKmers.find(*kmer)->second : createVertex(*kmer);
        addEdge(prevVertex, nextVertex, std::shared_ptr<MultiSampleEdge>(new MultiSampleEdge(isRef, count, numPruningSamples)));
        delete kmer;
        return nextVertex;
    }

}

void ReadThreadingGraph::threadSequence(SequenceForKmers & sequenceForKmers) {
    int uniqueStartPos = findStart(sequenceForKmers);
    if(uniqueStartPos == -1)
        return;

    std::shared_ptr<MultiDeBruijnVertex> startingVertex = getOrCreateKmerVertex(sequenceForKmers.sequence, uniqueStartPos);

    if(INCREASE_COUNTS_BACKWARDS) {
        increaseCountsInMatchedKmers(sequenceForKmers, startingVertex, startingVertex->getSequence(), kmerSize - 2);
    }

    if(sequenceForKmers.isRef) {
        if(refSource.getBases() != nullptr)
            throw std::invalid_argument("Found two refSources! prev:");

        refSource = Kmer(sequenceForKmers.sequence, sequenceForKmers.start, kmerSize);
    }

    std::shared_ptr<MultiDeBruijnVertex> vertex = startingVertex;
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

std::shared_ptr<MultiDeBruijnVertex> ReadThreadingGraph::getOrCreateKmerVertex(std::shared_ptr<uint8_t[]>sequence, const int start) {
    Kmer kmer(std::move(sequence), start, kmerSize);
    return getUniqueKmerVertex(kmer, true) ? uniqueKmers.find(kmer)->second : createVertex(kmer);
}

void ReadThreadingGraph::increaseCountsInMatchedKmers(SequenceForKmers & seqForKmers, const std::shared_ptr<MultiDeBruijnVertex>& vertex,
                                                      const std::shared_ptr<uint8_t[]>& originalKmer, int offset) {
    if(offset == -1)
        return;

    std::unordered_set<std::shared_ptr<MultiSampleEdge>> incomingEdges = incomingEdgesOf(vertex);
    for(std::unordered_set<std::shared_ptr<MultiSampleEdge>>::iterator iter = incomingEdges.begin(); iter != incomingEdges.end(); iter++) {
        std::shared_ptr<MultiDeBruijnVertex> prev = getEdgeSource(*iter);
        uint8_t suffix = prev->getSuffix();
        uint8_t seqBase = originalKmer.get()[offset];
        if(suffix == seqBase && (increaseCountsThroughBranches || inDegreeOf(vertex) == 1)) {
            (*iter)->incMultiplicity(seqForKmers.count);
            increaseCountsInMatchedKmers(seqForKmers, prev, originalKmer, offset-1);
        }
    }
}

void ReadThreadingGraph::buildGraphIfNecessary() {
    if(alreadyBuilt)
        return;

    //test
//    for(std::pair<std::string, std::vector<SequenceForKmers>> iter : pending) {
//        for(SequenceForKmers kmer : iter.second) {
//            uint8_t * sequence = kmer.sequence.get();
//            for(int i = 0; i < kmer.stop - kmer.start; i++) {
//                std::cout << sequence[i];
//            }
//            std::cout << std::endl;
//        }
//    }

    nonUniqueKmers = determineKmerSizeAndNonUniques(kmerSize, kmerSize);


    for(std::map<std::string, std::vector<SequenceForKmers>>::iterator miter = pending.begin(); miter != pending.end(); miter++){
        for(std::vector<SequenceForKmers>::iterator viter = miter->second.begin(); viter != miter->second.end(); viter++) {
            threadSequence(*viter);
        }
        std::unordered_map<std::shared_ptr<MultiSampleEdge>, IntrusiveEdge<MultiDeBruijnVertex>>::iterator eiter;
        for(eiter = edgeMap.begin(); eiter != edgeMap.end(); eiter++) {
            (*eiter->first).flushSingleSampleMultiplicity();
        }
    }

//    std::map<std::string, std::vector<SequenceForKmers>>::iterator miter = pending.begin();
//
//    for(std::vector<SequenceForKmers>::iterator viter = miter->second.begin(); viter != miter->second.end(); viter++) {
//        threadSequence(*viter);
//        std::cout << getVertexSet().size() << ", "<< getEdgeSet().size() << std::endl;
//    }
//    std::unordered_map<std::shared_ptr<MultiSampleEdge>, IntrusiveEdge<MultiDeBruijnVertex>>::iterator eiter;
//    for(eiter = edgeMap.begin(); eiter != edgeMap.end(); eiter++) {
//        (*eiter->first).flushSingleSampleMultiplicity();
//    }
//
//    miter++;
//    miter++;
//
//    if(miter != pending.end()) {
//        for(std::vector<SequenceForKmers>::iterator viter = miter->second.begin(); viter != miter->second.end(); viter++) {
//            threadSequence(*viter);
//            std::cout << getVertexSet().size() << ", "<< getEdgeSet().size() << std::endl;
//        }
//        for(eiter = edgeMap.begin(); eiter != edgeMap.end(); eiter++) {
//            (*eiter->first).flushSingleSampleMultiplicity();
//        }
//    }
//
//
//    miter--;
//
//    for(std::vector<SequenceForKmers>::iterator viter = miter->second.begin(); viter != miter->second.end(); viter++) {
//        threadSequence(*viter);
//        std::cout << getVertexSet().size() << ", "<< getEdgeSet().size() << std::endl;
//    }
//    for(eiter = edgeMap.begin(); eiter != edgeMap.end(); eiter++) {
//        (*eiter->first).flushSingleSampleMultiplicity();
//    }


    pending.clear();
    alreadyBuilt = true;

    //test
//    for(std::shared_ptr<MultiDeBruijnVertex> multiDeBruijnVertex : DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::getVertexSet()) {
//        uint8_t * sequence = multiDeBruijnVertex->getSequence().get();
//        for(int i = 0; i < multiDeBruijnVertex->getLength(); i++) {
//            std::cout << sequence[i];
//        }
//        std::cout << std::endl;
//    }

    for(std::map<Kmer, std::shared_ptr<MultiDeBruijnVertex>>::iterator kiter = uniqueKmers.begin(); kiter != uniqueKmers.end(); kiter++) {
        kiter->second->setAdditionalInfo(kiter->second->getAdditionalInfo() + '+');
    }

}

void ReadThreadingGraph::setThreadingStartOnlyAtExistingVertex(bool value) {
    startThreadingOnlyAtExistingVertex = value;
}

//void ReadThreadingGraph::setPending() {
//    pending.clear();
//    uint8_t* byte997137 = new uint8_t[99]{67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84};
//    SequenceForKmers sequenceForKmers7135 = {.name = "HWI-ST729_110151799:2:41:9503:140222_0_99", .sequence = byte997137, .start = 0, .stop = 99, .count = 1, .isRef = false};
//    std::vector<SequenceForKmers> v3;
//    v3.emplace_back(sequenceForKmers7135);
//    std::string key1491293244 = "H_LS-E2-A15C-01A-31D-A12B-09";
//
//    uint8_t* byte1007126 = new uint8_t[100]{67, 84, 65, 65, 67, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67};
//    SequenceForKmers sequenceForKmers7123 = {.name = "HWI-ST729_110151799:3:66:3953:197183_0_100", .sequence = byte1007126, .start = 0, .stop = 100, .count = 1, .isRef = false};
//    uint8_t* byte907131 = new uint8_t[90]{65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 67, 67, 65};
//    SequenceForKmers sequenceForKmers7124 = {.name = "HWI-ST729_110151799:3:64:7208:158288_0_90", .sequence = byte907131, .start = 0, .stop = 90, .count = 1, .isRef = false};
//    std::vector<SequenceForKmers> v2;
//    v2.emplace_back(sequenceForKmers7123);
//    v2.emplace_back(sequenceForKmers7124);
//    std::string key1617703971 = "H_LS-E2-A15C-10A-01D-A12B-09";
//
//    uint8_t* byte2867095 = new uint8_t[286]{67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 65, 67, 67, 67, 84, 65, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67, 67, 84, 65, 65, 67, 67};
//    SequenceForKmers sequenceForKmers7094 = {.name = "ref", .sequence = byte2867095, .start = 0, .stop = 286, .count = 1, .isRef = true};
//    std::vector<SequenceForKmers> v1;
//    v1.emplace_back(sequenceForKmers7094);
//    std::string key769618002 = "XXX_UNNAMED_XXX";
//    pending.insert(std::pair<std::string, std::vector<SequenceForKmers>>(key1617703971, v2));
//    pending.insert(std::pair<std::string, std::vector<SequenceForKmers>>(key769618002, v3));
//    pending.insert(std::pair<std::string, std::vector<SequenceForKmers>>(key1491293244, v1));
//}

bool ReadThreadingGraph::removeVertex(std::shared_ptr<MultiDeBruijnVertex> V) {
    std::shared_ptr<uint8_t[]> sequence(new uint8_t[V->getLength()]);
    memcpy(sequence.get(), V->getSequence().get(), V->getLength());
    bool result = DirectedSpecifics::removeVertex(V);

    if(result) {
        Kmer kmer(sequence, 0 ,kmerSize);
        uniqueKmers.erase(kmer);
    }

    return result;
}

void ReadThreadingGraph::removeSingletonOrphanVertices() {
    std::vector<std::shared_ptr<MultiDeBruijnVertex>> toRemove;
    std::unordered_set<std::shared_ptr<MultiDeBruijnVertex>> allvertex = getVertexSet();
    typename std::unordered_set<std::shared_ptr<MultiDeBruijnVertex>>::iterator viter;
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
    for ( std::shared_ptr<MultiDeBruijnVertex> v :  getVertexSet()) {
        if(outDegreeOf(v) == 0 && !isRefSink(v)) {
            attempted++;
            nRecovered += recoverDanglingTail(v, pruneFactor, minDanglingBranchLength, recoverAll);
        }
    }
}

int ReadThreadingGraph::recoverDanglingTail(std::shared_ptr<MultiDeBruijnVertex> vertex, int pruneFactor, int minDanglingBranchLength,
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
ReadThreadingGraph::generateCigarAgainstDownwardsReferencePath(std::shared_ptr<MultiDeBruijnVertex> vertex, int pruneFactor,
                                                               int minDanglingBranchLength, bool recoverAll) {
    int minTailPathLength = std::max(1, minDanglingBranchLength);
    std::deque<std::shared_ptr<MultiDeBruijnVertex>> altPath = findPathUpwardsToLowestCommonAncestor(vertex, pruneFactor, !recoverAll);
    if(altPath.empty() || isRefSource(altPath.front()) || altPath.size() < minTailPathLength + 1) {
        return nullptr;
    }
    std::shared_ptr<MultiDeBruijnVertex> toBlacklistedEdge =  altPath.size() > 1 ? altPath.at(1) : nullptr;
    std::vector<std::shared_ptr<MultiDeBruijnVertex>> refPath = getReferencePath(altPath.at(0), downwards, getHeaviestIncomingEdge(toBlacklistedEdge));
    std::vector<std::shared_ptr<MultiDeBruijnVertex>> newaltPath;
    for(std::shared_ptr<MultiDeBruijnVertex> vertex1 : altPath) {
        newaltPath.emplace_back(vertex1);
    }
    int refLength;
    std::shared_ptr<uint8_t[]> refBases = getBasesForPath(refPath, refLength, false);
    int altLength;
    std::shared_ptr<uint8_t[]> altBases = getBasesForPath(newaltPath, altLength, false);
    SWNativeAlignerWrapper wrapper = SWNativeAlignerWrapper();
    SmithWatermanAlignment* alignment = wrapper.align(refBases, refLength, altBases, altLength, &SmithWatermanAligner::STANDARD_NGS, LEADING_INDEL);
    return new DanglingChainMergeHelper(newaltPath, refPath, altBases, altLength, refBases, refLength, AlignmentUtils::removeTrailingDeletions(alignment->getCigar()));
}



std::deque<std::shared_ptr<MultiDeBruijnVertex>>
ReadThreadingGraph::findPathUpwardsToLowestCommonAncestor(std::shared_ptr<MultiDeBruijnVertex> vertex, int pruneFactor,
                                                          bool giveUpAtBranch) {
    std::deque<std::shared_ptr<MultiDeBruijnVertex>> ret;
    std::shared_ptr<MultiDeBruijnVertex> v = std::move(vertex);
    if(giveUpAtBranch) {
        while(!(inDegreeOf(v) != 1 || outDegreeOf(v) >= 2)) {
            std::shared_ptr<MultiSampleEdge> edge = incomingEdgeOf(v);
            if(edge->getPruningMultiplicity() < pruneFactor) {
                ret.clear();
            }
            else {
                ret.push_front(v);
            }
            v = getEdgeSource(edge);
        }
        ret.push_front(v);

        return outDegreeOf(v) > 1 ? ret : std::deque<std::shared_ptr<MultiDeBruijnVertex>>();
    } else {
        while(!(hasIncidentRefEdge(v) || inDegreeOf(v) == 0)) {
            std::shared_ptr<MultiSampleEdge> edge = getHeaviestIncomingEdge(v);
            if(edge->getPruningMultiplicity() < pruneFactor) {
                ret.clear();
            }
            else {
                ret.push_front(v);
            }
            v = getEdgeSource(edge);
        }
        ret.push_front(v);

        return outDegreeOf(v) > 1 && hasIncidentRefEdge(v) ? ret : std::deque<std::shared_ptr<MultiDeBruijnVertex>>();
    }
}

bool ReadThreadingGraph::hasIncidentRefEdge(std::shared_ptr<MultiDeBruijnVertex> v) {
    for(const std::shared_ptr<MultiSampleEdge>& edge : incomingEdgesOf(v)) {
        if(edge->getIsRef()) {
            return true;
        }
    }
    return false;
}

std::shared_ptr<MultiSampleEdge> ReadThreadingGraph::getHeaviestIncomingEdge(std::shared_ptr<MultiDeBruijnVertex> v) {
    std::unordered_set<std::shared_ptr<MultiSampleEdge>> incomingEdges = incomingEdgesOf(v);
    std::shared_ptr<MultiSampleEdge> ret;
    ret = *incomingEdges.begin();
    for(std::unordered_set<std::shared_ptr<MultiSampleEdge>>::iterator iter = incomingEdges.begin(); iter != incomingEdges.end(); iter++) {
        if(ret->getPruningMultiplicity() < (*iter)->getPruningMultiplicity())
            ret = *iter;
    }
    return ret;
}

std::vector<std::shared_ptr<MultiDeBruijnVertex>>
ReadThreadingGraph::getReferencePath(std::shared_ptr<MultiDeBruijnVertex>start, TraversalDirection direction,
                                     const std::shared_ptr<MultiSampleEdge>& blacklistedEdge) {
    std::vector<std::shared_ptr<MultiDeBruijnVertex>> path;

    std::shared_ptr<MultiDeBruijnVertex> v = std::move(start);

    while (v != nullptr) {
        path.emplace_back(v);
        v = (direction == downwards ? getNextReferenceVertex(v, true, blacklistedEdge) : getPrevReferenceVertex(v));
    }
    return path;
}

std::shared_ptr<uint8_t[]>ReadThreadingGraph::getBasesForPath(const std::vector<std::shared_ptr<MultiDeBruijnVertex>>& path, int &length, bool expandSource) {
    int tmpLength = 300;
    int start = 0;
    std::shared_ptr<uint8_t[]> tmp(new uint8_t[300]);
    for(const std::shared_ptr<MultiDeBruijnVertex>& v : path) {
        if(expandSource && isSource(v)) {
            std::shared_ptr<uint8_t[]> seq = v->getSequence();
            int seqLength = v->getLength();
            while(start + seqLength > tmpLength) {
                tmpLength *= 2;
                std::shared_ptr<uint8_t[]> newtmp(new uint8_t[tmpLength]);
                memcpy(newtmp.get(), tmp.get(), tmpLength);
                tmp = newtmp;
            }
            for(int i = 0; i < seqLength; i++) {
                tmp.get()[start + i] = seq.get()[seqLength-1-i];
            }
            start += seqLength;
        } else {
            while(start + 1 > tmpLength) {
                tmpLength *= 2;
                std::shared_ptr<uint8_t[]> newtmp(new uint8_t[tmpLength]);
                memcpy(newtmp.get(), tmp.get(), tmpLength);
                tmp = newtmp;
            }
            tmp.get()[start] = v->getSuffix();
            start++;
        }
    }
    std::shared_ptr<uint8_t[]> ret(new uint8_t[start+1]);
    memcpy(ret.get(), tmp.get(), start);
    ret.get()[start] = '\0';
    length = start;
    return ret;
}

bool ReadThreadingGraph::cigarIsOkayToMerge(std::shared_ptr<Cigar> &cigar, bool requireFirstElementM, bool requireLastElementM) {
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
    addEdge(danglingTailMergeResult->danglingPath.at(altIndexToMerge), danglingTailMergeResult->referencePath.at(refIndexToMerge), std::shared_ptr<MultiSampleEdge>(new MultiSampleEdge(
            false, 1, numPruningSamples)));
    return 1;
}

int ReadThreadingGraph::longestSuffixMatch(std::shared_ptr<uint8_t[]>seq, int seqLength, std::shared_ptr<uint8_t[]> kmer, int kmerLength, int seqStart) {
    for (int len = 1; len <= kmerLength; len++) {
        int seqI = seqStart - len + 1;
        int kmerI = kmerLength - len;
        if ( seqI < 0 || seq.get()[seqI] != kmer.get()[kmerI] ) {
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

    std::vector<std::shared_ptr<MultiDeBruijnVertex>> danglingHeads;
    for ( std::shared_ptr<MultiDeBruijnVertex> v :  getVertexSet()) {
        if( DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::inDegreeOf(v) == 0 && !isRefSource(v))
            danglingHeads.emplace_back(v);
    }

    int attempted = 0;
    int nRecovered = 0;
    for ( std::shared_ptr<MultiDeBruijnVertex> v :  danglingHeads) {
        if(outDegreeOf(v) == 0 && !isRefSink(v)) {
            attempted++;
            nRecovered += recoverDanglingHead(v, pruneFactor, minDanglingBranchLength, recoverAll);
        }
    }
}

int ReadThreadingGraph::recoverDanglingHead(std::shared_ptr<MultiDeBruijnVertex>vertex, int pruneFactor, int minDanglingBranchLength,
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
ReadThreadingGraph::generateCigarAgainstUpwardsReferencePath(std::shared_ptr<MultiDeBruijnVertex> vertex, int pruneFactor,
                                                             int minDanglingBranchLength, bool recoverAll) {
    std::deque<std::shared_ptr<MultiDeBruijnVertex>> altPath = findPathDownwardsToHighestCommonDescendantOfReference(vertex, pruneFactor, !recoverAll);
    if(altPath.empty() || isRefSink(altPath.front()) || altPath.size() < minDanglingBranchLength + 1) {
        return nullptr;
    }
    std::vector<std::shared_ptr<MultiDeBruijnVertex>> refPath = getReferencePath(altPath.at(0), upwards, nullptr);
    std::vector<std::shared_ptr<MultiDeBruijnVertex>> newaltPath;
    for(std::shared_ptr<MultiDeBruijnVertex> vertex1 : altPath) {
        newaltPath.emplace_back(vertex1);
    }
    int refLength;
    std::shared_ptr<uint8_t[]> refBases = getBasesForPath(refPath, refLength, false);
    int altLength;
    std::shared_ptr<uint8_t[]> altBases = getBasesForPath(newaltPath, altLength, false);
    SWNativeAlignerWrapper wrapper = SWNativeAlignerWrapper();
    SmithWatermanAlignment* alignment = wrapper.align(refBases, refLength, altBases, altLength, &SmithWatermanAligner::STANDARD_NGS, LEADING_INDEL);
    return new DanglingChainMergeHelper(newaltPath, refPath, altBases, altLength, refBases, refLength, AlignmentUtils::removeTrailingDeletions(alignment->getCigar()));
}

std::deque<std::shared_ptr<MultiDeBruijnVertex>>
ReadThreadingGraph::findPathDownwardsToHighestCommonDescendantOfReference(std::shared_ptr<MultiDeBruijnVertex>vertex, int pruneFactor,
                                                                          bool giveUpAtBranch) {
    std::deque<std::shared_ptr<MultiDeBruijnVertex>> ret;
    std::shared_ptr<MultiDeBruijnVertex> v = vertex;
    if(giveUpAtBranch) {
        while(isReferenceNode(v) || outDegreeOf(v) != 1) {
            std::shared_ptr<MultiSampleEdge> edge = outgoingEdgeOf(v);
            if(edge->getPruningMultiplicity() < pruneFactor) {
                ret.clear();
            }
            else {
                ret.push_front(v);
            }
            v = getEdgeTarget(edge);
        }
        ret.push_front(v);

        return isReferenceNode(v)? ret : std::deque<std::shared_ptr<MultiDeBruijnVertex>>();
    } else {
        while(isReferenceNode(v) || outDegreeOf(v) == 0) {
            std::shared_ptr<MultiSampleEdge> edge = getHeaviestOutgoingEdge(v);
            if(edge->getPruningMultiplicity() < pruneFactor) {
                ret.clear();
            }
            else {
                ret.push_front(v);
            }
            v = getEdgeTarget(edge);
        }
        ret.push_front(v);

        return isReferenceNode(v) ? ret : std::deque<std::shared_ptr<MultiDeBruijnVertex>>();
    }
}

std::shared_ptr<MultiSampleEdge> ReadThreadingGraph::getHeaviestOutgoingEdge(std::shared_ptr<MultiDeBruijnVertex> v) {
    std::unordered_set<std::shared_ptr<MultiSampleEdge>> outgoing = outgoingEdgesOf(v);
    std::shared_ptr<MultiSampleEdge> ret;
    ret = *outgoing.begin();
    for(std::unordered_set<std::shared_ptr<MultiSampleEdge>>::iterator iter = outgoing.begin(); iter != outgoing.end(); iter++) {
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

    addEdge(danglingHeadMergeResult->referencePath.at(indexesToMerge+1), danglingHeadMergeResult->danglingPath.at(indexesToMerge),
            std::shared_ptr<MultiSampleEdge>(new MultiSampleEdge(
                    false, 1, numPruningSamples)) );

    return 1;
}

int ReadThreadingGraph::bestPrefixMatch(const std::shared_ptr<uint8_t[]>&path1, const std::shared_ptr<uint8_t[]>path2, int maxIndex) {
    int maxMismatches = getMaxMismatches(maxIndex);
    int mismatches = 0;
    int index = 0;
    int lastGoodIndex = -1;
    uint8_t * path1_ = path1.get();
    uint8_t * path2_ = path2.get();
    while ( index < maxIndex ) {
        if ( path1_[index] != path2_[index] ) {
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

    std::shared_ptr<MultiDeBruijnVertex> danglingSource = danglingHeadMergeResult->danglingPath.at(indexOfLastDanglingNode);
    danglingHeadMergeResult->danglingPath.erase(danglingHeadMergeResult->danglingPath.begin() + indexOfLastDanglingNode);
    std::shared_ptr<uint8_t[]> danglingSourceSeq = danglingSource->getSequence();
    int danglingSourceSeqLength = danglingSource->getLength();
    std::shared_ptr<uint8_t[]> refSourceSequence = danglingHeadMergeResult->referencePath.at(indexOfRefNodeToUse)->getSequence();
    std::shared_ptr<uint8_t[]> sequenceToExtend(new uint8_t[numNodesToExtend + danglingSourceSeqLength]);
    uint8_t * sequenceToExtend_ = sequenceToExtend.get();
    uint8_t * refSourceSequence_ = refSourceSequence.get();
    for ( int i = 0; i < numNodesToExtend; i++ ) {
        sequenceToExtend_[i] = refSourceSequence_[i];
    }
    memcpy(sequenceToExtend.get()+numNodesToExtend, danglingSourceSeq.get(), danglingSourceSeqLength);

    std::shared_ptr<MultiSampleEdge> sourceEdge = getHeaviestOutgoingEdge(danglingSource);
    std::shared_ptr<MultiDeBruijnVertex> prevV = getEdgeTarget(sourceEdge);
    std::shared_ptr<MultiSampleEdge> ret = removeEdge(danglingSource, prevV);
    for( int i = numNodesToExtend; i > 0; i-- ) {
        std::shared_ptr<uint8_t[]> tmp(new uint8_t[kmerSize]);
        memcpy(tmp.get(), sequenceToExtend.get() + i, kmerSize);
        std::shared_ptr<MultiDeBruijnVertex> newV = std::shared_ptr<MultiDeBruijnVertex>(new MultiDeBruijnVertex(tmp, kmerSize));
        addVertex(newV);
        std::shared_ptr<MultiSampleEdge> newE = addEdge(newV, prevV);
        newE->setMultiplicity(sourceEdge->getMultiplicity());
        danglingHeadMergeResult->danglingPath.emplace_back(newV);
        prevV = newV;
    }
    return true;
}

std::shared_ptr<MultiSampleEdge> ReadThreadingGraph::createEdge(std::shared_ptr<MultiDeBruijnVertex>, std::shared_ptr<MultiDeBruijnVertex> ) {
    return std::shared_ptr<MultiSampleEdge>(new MultiSampleEdge(false, 1, numPruningSamples));
}

std::shared_ptr<SeqGraph> ReadThreadingGraph::toSequenceGraph() {
    buildGraphIfNecessary();
    std::shared_ptr<SeqGraph> seqGraph(new SeqGraph(kmerSize));
    std::map<std::shared_ptr<MultiDeBruijnVertex>, std::shared_ptr<SeqVertex>> vertexMap;
    for(std::shared_ptr<MultiDeBruijnVertex> dv : DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::getVertexSet()) {
        std::shared_ptr<SeqVertex> sv(new SeqVertex(dv->getAdditionalSequence(
                DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::isSource(dv)), dv->getAdditionalLength(DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::isSource(dv))));
        sv->setAdditionalInfo(dv->getAdditionalInfo());
        vertexMap.insert(std::pair<std::shared_ptr<MultiDeBruijnVertex>, std::shared_ptr<SeqVertex>>(dv, sv));
        seqGraph->addVertex(sv);
    }

    std::unordered_map<std::shared_ptr<MultiSampleEdge>, IntrusiveEdge<MultiDeBruijnVertex>>::iterator eiter;
    for(eiter = edgeMap.begin(); eiter != edgeMap.end(); eiter++) {
        std::shared_ptr<SeqVertex> seqInV = vertexMap.at(getEdgeSource(eiter->first));
        std::shared_ptr<SeqVertex> seqOutV = vertexMap.at(getEdgeTarget(eiter->first));
        seqGraph->addEdge(seqInV, seqOutV, std::shared_ptr<BaseEdge>(new BaseEdge(eiter->first->getIsRef(), eiter->first->getMultiplicity())));
    }
    return seqGraph;
}

ReadThreadingGraph::ReadThreadingGraph(int kmerSize, bool debugGraphTransformations,
                                       uint8_t minBaseQualityToUseInAssembly, int numPruningSamples) : kmerSize(kmerSize), minBaseQualityToUseInAssembly(minBaseQualityToUseInAssembly), debugGraphTransformations(debugGraphTransformations),
                                                                                                       refSource(Kmer(nullptr, 0)), numPruningSamples(numPruningSamples){
    Mutect2Utils::validateArg(kmerSize > 0, "bad minkKmerSize");
    resetToInitialState();
}

void ReadThreadingGraph::resetToInitialState() {
    pending.clear();
    nonUniqueKmers.clear();
    uniqueKmers.clear();
    alreadyBuilt = false;
}

void ReadThreadingGraph::addSequence(std::string seqName, std::shared_ptr<uint8_t[]>sequence, int length, int count, bool isRef) {
    addSequence(std::move(seqName), ANONYMOUS_SAMPLE, sequence, 0, length, count, isRef);
}

void ReadThreadingGraph::addSequence(std::string seqName, std::shared_ptr<uint8_t[]> sequence, int length, bool isRef) {
    addSequence(std::move(seqName), sequence, length, 1, isRef);
}

bool ReadThreadingGraph::isLowComplexity() {
    return nonUniqueKmers.size() * 4 > uniqueKmers.size();
}

bool ReadThreadingGraph::hasCycles() {
    DFS_CycleDetect<MultiDeBruijnVertex,MultiSampleEdge> detect = DFS_CycleDetect<MultiDeBruijnVertex,MultiSampleEdge>(*this);
    return detect.detectCycles();
}







