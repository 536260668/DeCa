//
// Created by cluster on 22-9-27.
//

#ifndef MUTECT2CPP_MASTER_TIRETREENODE_H
#define MUTECT2CPP_MASTER_TIRETREENODE_H

#include <vector>

class tireTreeNode {
private:
    int size;
    std::vector<int> index;
    std::vector<tireTreeNode *> childs;

public:
    tireTreeNode();
    ~tireTreeNode();
    tireTreeNode(const std::vector<int> &index);
    void addChild(tireTreeNode *node);
    std::vector<tireTreeNode *> getChild();
    tireTreeNode & operator=(tireTreeNode const& node);
    void addIndex(int child);
    std::vector<int> getIndex() const;
};


#endif //MUTECT2CPP_MASTER_TIRETREENODE_H
