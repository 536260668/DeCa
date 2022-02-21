//
// Created by 梦想家xixi on 2021/11/19.
//

#include "MergeDiamonds.h"
#include "SharedVertexSequenceSplitter.h"

bool MergeDiamonds::tryToTransform(std::shared_ptr<SeqVertex> top) {
    Mutect2Utils::validateArg(top != nullptr, "Null is not allowed there");
    std::shared_ptr<SeqGraph> graph1 = getGraph();
    std::unordered_set<std::shared_ptr<SeqVertex>> middles = graph1->outgoingVerticesOf(top);
    if(middles.size() <= 1) {
        return false;
    }

    std::shared_ptr<SeqVertex> bottom = nullptr;
    for(std::shared_ptr<SeqVertex> mi : middles) {
        if(graph1->outDegreeOf(mi) < 1) {
            return false;
        }
        if(graph1->inDegreeOf(mi) != 1) {
            return false;
        }

        for(std::shared_ptr<SeqVertex> mt : graph1->outgoingVerticesOf(mi)) {
            if(bottom == nullptr) {
                bottom = mt;
            } else if (bottom != mt) {
                return false;
            }
        }
    }

    if(graph1->inDegreeOf(bottom) != middles.size()) {
        return false;
    }
    if(getDontModifyGraphEvenIfPossible()) {
        return true;
    }
    SharedVertexSequenceSplitter splitter = SharedVertexSequenceSplitter(getGraph(), middles);
    return splitter.meetsMinMergableSequenceForEitherPrefixOrSuffix(1) && splitter.splitAndUpdate(top, bottom);
}
