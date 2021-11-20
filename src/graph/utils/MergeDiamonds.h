//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_MERGEDIAMONDS_H
#define MUTECT2CPP_MASTER_MERGEDIAMONDS_H

#include "VertexBasedTransformer.h"

class MergeDiamonds : public VertexBasedTransformer{
public:
    MergeDiamonds(SeqGraph* seqGraph) : VertexBasedTransformer(seqGraph) {}

protected:
    bool tryToTransform(SeqVertex* top) override;
};


#endif //MUTECT2CPP_MASTER_MERGEDIAMONDS_H
