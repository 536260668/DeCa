//
// Created by 梦想家xixi on 2021/11/24.
//

#include "KBestHaplotype.h"
#include <cmath>

KBestHaplotype::KBestHaplotype(std::shared_ptr<SeqVertex> initialVertex, std::shared_ptr<DirectedSpecifics<SeqVertex, BaseEdge>> graph) : Path<SeqVertex, BaseEdge>(initialVertex, graph){
    score = 0;
}

KBestHaplotype::KBestHaplotype(std::shared_ptr<KBestHaplotype> p, std::shared_ptr<BaseEdge> edge, int totalOutgoingMultiplicity) : Path<SeqVertex, BaseEdge>(*p, edge){
    score = p->getScore() + std::log10(edge->getMultiplicity()) - std::log10(totalOutgoingMultiplicity);
    isReference &= edge->getIsRef();
}

std::shared_ptr<Haplotype> KBestHaplotype::getHaplotype() {
    int length = 0;
    std::shared_ptr<uint8_t[]> base = getBases(length);
    std::shared_ptr<Haplotype> haplotype(new Haplotype(base, length, getIsReference()));
    haplotype->setScore(score);
    return haplotype;
}
