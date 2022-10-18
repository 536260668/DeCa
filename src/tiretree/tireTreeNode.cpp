//
// Created by cluster on 22-9-27.
//

#include <cstring>
#include "tireTreeNode.h"
#include <iostream>

tireTreeNode::tireTreeNode() {
    size = 0;
}

tireTreeNode &tireTreeNode::operator=(tireTreeNode const&node) {
    size = node.size;
    index = node.index;
    childs = node.childs;
    return *this;
}


std::vector<int>& tireTreeNode::getIndex(){
    return index;
}

void tireTreeNode::addIndex(int child) {
    index.emplace_back(child);
}

tireTreeNode::~tireTreeNode() {

}

tireTreeNode::tireTreeNode(const std::vector<int> &index) {
    this->index = index;
    size = index.size();
}

void tireTreeNode::addChild(tireTreeNode *node) {
    childs.emplace_back(node);
}

std::vector<tireTreeNode *>& tireTreeNode::getChild() {
    return childs;
}



