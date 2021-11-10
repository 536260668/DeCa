//
// Created by 梦想家xixi on 2021/10/18.
//

#include "ReadThreadingGraph.h"
#include <string>
#include <utility>
#include "Mutect2Utils.h"

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



