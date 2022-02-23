//
// Created by 梦想家xixi on 2021/10/20.
//

#ifndef MUTECT2CPP_MASTER_MULTIDEBRUIJNVERTEX_H
#define MUTECT2CPP_MASTER_MULTIDEBRUIJNVERTEX_H

#include "BaseVertex.h"
#include <vector>

class MultiDeBruijnVertex : public BaseVertex{
private:
    std::vector<std::string> read;
    long hashCode;
    bool mergeIdenticalNodes;

public:
    /**
    * Create a new MultiDeBruijnVertex with kmer sequence
    * @param mergeIdenticalNodes should nodes with the same sequence be treated as equal?
    * @param sequence the kmer sequence
    */
    MultiDeBruijnVertex(std::shared_ptr<uint8_t[]> sequence, int length, bool mergeIdenticalNodes);

    MultiDeBruijnVertex(std::shared_ptr<uint8_t[]> sequence, int length);

    ~MultiDeBruijnVertex() override = default;

    bool operator==(const MultiDeBruijnVertex &other) const;

    bool operator<(const MultiDeBruijnVertex &other) const;

    int getKmerSize() const {return getLength();}

    uint8_t getSuffix() const {
        return  sequence.get()[getKmerSize() - 1];
    }

    std::shared_ptr<uint8_t[]> getAdditionalSequence(bool source) override;

    std::shared_ptr<uint8_t[]> getSuffixAsArray() const;

    int getAdditionalLength(bool source);

    int getAdditionalSequenceLength(bool isSource) override;
};


#endif //MUTECT2CPP_MASTER_MULTIDEBRUIJNVERTEX_H
