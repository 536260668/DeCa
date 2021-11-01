//
// Created by 梦想家xixi on 2021/10/18.
//

#include "ReadThreadingAssembler.h"
#include <string>
#include <utility>
#include "Mutect2Utils.h"

void ReadThreadingAssembler::addRead(SAMRecord read) {
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

bool ReadThreadingAssembler::baseIsUsableForAssembly(uint8_t base, uint8_t qual) const {
    return base != 'N' && qual >= minBaseQualityToUseInAssembly;
}

void
ReadThreadingAssembler::addSequence(std::string seqName, std::string& sampleName, const uint8_t *sequence, int start, int stop,
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

std::vector<Kmer> ReadThreadingAssembler::determineNonUniqueKmers(SequenceForKmers &sequenceForKmers, const int kmerSize) {
    std::set<Kmer> allKmers;
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

std::vector<SequenceForKmers> ReadThreadingAssembler::getAllPendingSequences(){
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

std::set<Kmer> ReadThreadingAssembler::determineKmerSizeAndNonUniques(const int minKmerSize, const int maxKmerSize) {
    std::vector<SequenceForKmers> withNonUniques = getAllPendingSequences();
    std::set<Kmer> nonUniqueKmers_m;

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

MultiDeBruijnVertex ReadThreadingAssembler::createVertex(Kmer kmer) {
    MultiDeBruijnVertex newVertex(kmer.getBases(), kmer.getLength(), false);
    unsigned prevSize = unmodifiableVertexSet.size();
    if(!unmodifiableVertexSet.insert(newVertex).second)
        throw std::invalid_argument("Adding vertex to graph didn't increase the graph size");

    if(nonUniqueKmers.find(kmer) == nonUniqueKmers.end() && uniqueKmers.find(kmer) == uniqueKmers.end())
        uniqueKmers.insert(std::pair<Kmer, MultiDeBruijnVertex>(kmer, newVertex));

    return newVertex;
}

void ReadThreadingAssembler::test(Kmer kmer) {
    createVertex(kmer);
}
