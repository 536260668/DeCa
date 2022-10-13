//
// Created by cluster on 22-9-27.
//

#include "buildTreeUtils.h"
#include "intel/common/avx.h"
#include <iostream>
#include <deque>

tireTreeNode* buildTreeUtils::buildTreeWithHaplotype(const std::vector<std::shared_ptr<Haplotype>> &haplotypes, bool isFloat) {
    if(haplotypes.empty()) {
        return nullptr;
    }
    char *reference = nullptr;
    int referenceLength = 0;
    int record;
    for(int j = 0; j < haplotypes.size(); j++) {
        if(haplotypes[j]->getIsReference()) {
            reference = reinterpret_cast<char*>(haplotypes[j]->getBases().get());
            referenceLength = haplotypes[j]->getBasesLength();
            record = j;
            break;
        }
    }
    if(reference == nullptr) {
        throw std::invalid_argument("there is no refHaplotypes");
    }
    int avxLen = isFloat ? sizeof(float)*8 : sizeof(double)*8;
    int i = 0;
    tireTreeNode *root = new tireTreeNode();
    tireTreeNode *father = root;
    while(referenceLength > i * avxLen) {
        i++;
        tireTreeNode *child = nullptr;
        child = new tireTreeNode({record});
        father->addChild(child);
        father = child;
    }
    for(int j = 0; j < haplotypes.size(); j++) {
        if(haplotypes[j]->getIsReference()) {
            continue;
        }
        int baseLen = haplotypes[j]->getBasesLength();
        char* bases = reinterpret_cast<char*>(haplotypes[j]->getBases().get());
        i = 0;
        father = root;
        while(i * avxLen < baseLen) {
            i++;
            bool flag = false;
            tireTreeNode *tmp;
            for(const auto & node : father->getChild()) {
                if(i * avxLen < baseLen) {
                    char *nodebases = reinterpret_cast<char*>(haplotypes[node->getIndex()[0]]->getBases().get());
                    int nodebaseLen =  haplotypes[node->getIndex()[0]]->getBasesLength();
                    if(i * avxLen < nodebaseLen && isEqual(nodebases+(i-1)*avxLen, bases+(i-1)*avxLen, avxLen)) {
                        flag = true;
                        tmp = node;
                        break;
                    }
                } else {
                    char *nodebases = reinterpret_cast<char*>(haplotypes[node->getIndex()[0]]->getBases().get());
                    int nodebaseLen =  haplotypes[node->getIndex()[0]]->getBasesLength();
                    if(baseLen == nodebaseLen && isEqual(nodebases+(i-1)*avxLen, bases+(i-1)*avxLen, baseLen-avxLen*(i-1))) {
                        throw std::invalid_argument("there are two same haplotypes");
                    }
                }
            }
            if(!flag) {
                tireTreeNode *child = new tireTreeNode({j});
                father->addChild(child);
                father = child;
            } else {
                tmp->addIndex(j);
                father = tmp;
            }
        }
    }
    return root;
}

size_t buildTreeUtils::avxLength() {
    return 32;
}

bool buildTreeUtils::isEqual(char *c1, char *c2, int len) {
    for(int i = 0; i < len; i++) {
        if(c1[i] != c2[i]) {
            return false;
        }
    }
    return true;
}

void buildTreeUtils::deleteTree(tireTreeNode *root) {
    if(root->getChild().empty()) {
        delete root;
    } else {
        for(auto & node : root->getChild()) {
            deleteTree(node);
        }
        delete root;
    }
}

void buildTreeUtils::printLayerTree(tireTreeNode *root) {
    std::deque<tireTreeNode *> records;
    records.emplace_back(root);
    while(!records.empty()) {
        std::deque<tireTreeNode *> tmp;
        while(!records.empty()) {
            tireTreeNode *node = records.front();
            records.pop_front();
            for(auto newnode : node->getChild()) {
                tmp.push_back(newnode);
                for(int i : newnode->getIndex()) {
                    std::cout << i << ", ";
                }
                std::cout << "    ";
            }
        }
        std::cout << std::endl;
        records = tmp;
    }
}
