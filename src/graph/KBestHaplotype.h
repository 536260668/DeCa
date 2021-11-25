//
// Created by 梦想家xixi on 2021/11/24.
//

#ifndef MUTECT2CPP_MASTER_KBESTHAPLOTYPE_H
#define MUTECT2CPP_MASTER_KBESTHAPLOTYPE_H

#include "path/Path.h"
#include "SeqGraph.h"
#include "haplotype/Haplotype.h"

class KBestHaplotype : public Path<SeqVertex, BaseEdge>{
private:
    double score;
    bool isReference;

public:
    double getScore() {return score;}
    bool getIsReference() {return isReference;}
    KBestHaplotype(SeqVertex* initialVertex, DirectedSpecifics<SeqVertex, BaseEdge> & graph);
    KBestHaplotype(KBestHaplotype* p, BaseEdge* edge, int totalOutgoingMultiplicity);
    Haplotype* getHaplotype();
};


#endif //MUTECT2CPP_MASTER_KBESTHAPLOTYPE_H
