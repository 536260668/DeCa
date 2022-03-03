//
// Created by 梦想家xixi on 2022/3/1.
//

#ifndef MUTECT2CPP_MASTER_HAPLOTYPEDATAHOLDER_H
#define MUTECT2CPP_MASTER_HAPLOTYPEDATAHOLDER_H

#include <cstdint>
#include <memory>
#include <utility>

struct HaplotypeDataHolder {
    std::shared_ptr<uint8_t[]> haplotypeBases;
    unsigned length;
    HaplotypeDataHolder(std::shared_ptr<uint8_t[]>  haplotypeBases, unsigned length) : haplotypeBases(std::move(haplotypeBases)), length(length) {};
};


#endif //MUTECT2CPP_MASTER_HAPLOTYPEDATAHOLDER_H
