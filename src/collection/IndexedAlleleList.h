//
// Created by 梦想家xixi on 2022/3/1.
//

#ifndef MUTECT2CPP_MASTER_INDEXEDALLELELIST_H
#define MUTECT2CPP_MASTER_INDEXEDALLELELIST_H

#include "IndexedSet.h"

template<typename E, typename _hash, typename _equal>
class IndexedAlleleList {
private:
    IndexedSet<E, _hash, _equal> alleles;

public:
    IndexedAlleleList() = default;

    /**
     * Constructs a new allele-list from a collection of alleles.
     *
     * <p>
     *     Repeats in the input collection will be ignored (keeping the first one). The order of alleles in the
     *     resulting list is the same as in the natural traversal of the input collection.
     *
     * </p>
     * @param alleles the original allele collection
     *
     * @throws IllegalArgumentException if {@code alleles} is {@code null} or contains {@code null}s.
     */
    IndexedAlleleList(const std::vector<E> & alleles) : alleles(alleles) {}

    int numberOfAlleles() const {
        return alleles.size();
    }

    int indexOfAllele(const E & allele) const {
        return alleles.indexOf(allele);
    }

    E getAllele(const int index) {
        return alleles.get(index);
    }
};

#endif //MUTECT2CPP_MASTER_INDEXEDALLELELIST_H
