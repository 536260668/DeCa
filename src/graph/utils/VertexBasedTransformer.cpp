//
// Created by 梦想家xixi on 2021/11/19.
//

#include "VertexBasedTransformer.h"

VertexBasedTransformer::VertexBasedTransformer(SeqGraph* graph) : graph(graph){
    Mutect2Utils::validateArg(graph != nullptr, "Null is not allowed there");
}

bool VertexBasedTransformer::transformUntilComplete() {
    bool didAtLeastOneTransform = false;
    bool foundNodesToMerge = true;
    while( foundNodesToMerge ) {
        foundNodesToMerge = false;
        for(const std::shared_ptr<SeqVertex>& v : graph->getSortedVertexList()) {
            foundNodesToMerge = tryToTransform(v);
            if(foundNodesToMerge) {
                didAtLeastOneTransform = true;
                break;
            }
        }

    }
    return didAtLeastOneTransform;
}
