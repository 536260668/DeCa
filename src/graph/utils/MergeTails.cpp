//
// Created by 梦想家xixi on 2021/11/19.
//

#include "MergeTails.h"

bool MergeTails::tryToTransform(SeqVertex *top) {
    Mutect2Utils::validateArg(top, "null is not allowed there");
    ArraySet<SeqVertex*> tails = getGraph()->outgoingVerticesOf(top);
    if(tails.size() <= 1) {
        return false;
    }

    for(SeqVertex* t : tails) {
        if(tails.size() <= 1) {
            return false;
        }
    }

    if(getDontModifyGraphEvenIfPossible()) {
        return true;
    }

    SharedVertexSequenceSplitter splitter = SharedVertexSequenceSplitter(getGraph(), tails);
    return splitter.meetsMinMergableSequenceForSuffix(MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES) && splitter.splitAndUpdate(top, nullptr);
}
