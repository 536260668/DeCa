//
// Created by 梦想家xixi on 2021/11/19.
//

#include "GraphUtils.h"
#include "Mutect2Utils.h"

int GraphUtils::minKmerLength(std::list<std::pair<uint8_t *, int>>&  kmers) {
    if(kmers.empty())
        return 0;
    int ret = kmers.begin()->second;
    for(std::pair<uint8_t*, int> kmer : kmers) {
        if(kmer.second < ret)
            ret = kmer.second;
    }
    return ret;
}

int GraphUtils::commonMaximumPrefixLength(std::list<std::pair<uint8_t *, int>> & listOfBytes) {
    int minLength = minKmerLength(listOfBytes);
    for (int i = 0; i < minLength; i++) {
        std::list<std::pair<uint8_t *, int>>::iterator liter = listOfBytes.begin();
        uint8_t b = listOfBytes.begin()->first[i];
        liter++;
        for(; liter != listOfBytes.end(); liter++) {
            if(b != liter->first[i]) {
                return i;
            }
        }
    }
    return minLength;
}

int GraphUtils::commonMaximumSuffixLength(std::list<std::pair<uint8_t *, int>> & listOfBytes,  const int minLength) {
    for(int suffixLen = 0; suffixLen < minLength; suffixLen++) {
        uint8_t b = listOfBytes.begin()->first[listOfBytes.begin()->second - suffixLen -1];
        std::list<std::pair<uint8_t *, int>>::iterator liter = listOfBytes.begin();
        liter++;
        for(; liter != listOfBytes.end(); liter++) {
            if(b != liter->first[liter->second - suffixLen - 1]) {
                return suffixLen;
            }
        }
    }
    return minLength;
}

std::list<std::pair<uint8_t *, int>> GraphUtils::getKmers(std::vector<SeqVertex *> vertices) {
    Mutect2Utils::validateArg(!vertices.empty(), "no vertex");
    std::list<std::pair<uint8_t *, int>> ret;
    for(SeqVertex* v : vertices) {
        ret.emplace_back(std::pair<uint8_t *, int>(v->getSequence(), v->getLength()));
    }
    return ret;
}
