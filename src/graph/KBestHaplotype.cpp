//
// Created by 梦想家xixi on 2021/11/24.
//

#include "KBestHaplotype.h"
#include <cmath>

KBestHaplotype::KBestHaplotype(SeqVertex *initialVertex, DirectedSpecifics<SeqVertex, BaseEdge> & graph) : Path<SeqVertex, BaseEdge>(initialVertex, graph){
    score = 0;
}

KBestHaplotype::KBestHaplotype(KBestHaplotype *p, BaseEdge *edge, int totalOutgoingMultiplicity) : Path<SeqVertex, BaseEdge>(*p, edge){
    score = p->getScore() + std::log10(edge->getMultiplicity()) - std::log10(totalOutgoingMultiplicity);
}

Haplotype *KBestHaplotype::getHaplotype() {
    int length = 0;
    uint8_t * base = getBases(length);
    Haplotype* haplotype = new Haplotype(base, length, getIsReference());
    haplotype->setScore(score);
    return haplotype;
}
