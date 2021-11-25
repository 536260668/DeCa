//
// Created by 梦想家xixi on 2021/11/20.
//

#include "MergeCommonSuffices.h"
#include "SharedSequenceMerger.h"
bool MergeCommonSuffices::tryToTransform(SeqVertex *bottom) {
    Mutect2Utils::validateArg(bottom, "null is not allowed there");
    return SharedSequenceMerger::merge(getGraph(), bottom);
}
