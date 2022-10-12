//
// Created by cluster on 22-9-27.
//

#ifndef MUTECT2CPP_MASTER_BUILDTREEUTILS_H
#define MUTECT2CPP_MASTER_BUILDTREEUTILS_H

#include "Haplotype.h"
#include "tireTreeNode.h"
#include <vector>

class buildTreeUtils {
public:
    static tireTreeNode* buildTreeWithHaplotype(const std::vector<std::shared_ptr<Haplotype>> &haplotypes, bool isFloat);
    static void printLayerTree(tireTreeNode *root);
    static void deleteTree(tireTreeNode *root);

private:
    static size_t avxLength();
    static bool isEqual(char *c1, char *c2, int len);
};


#endif //MUTECT2CPP_MASTER_BUILDTREEUTILS_H
