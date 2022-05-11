//
// Created by lhh on 5/7/22.
//

#ifndef MUTECT2CPP_MASTER_FRAGMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_FRAGMENTCOLLECTION_H

#include <vector>
#include <memory>
#include <cassert>
#include <string>
#include <unordered_map>
#include "samtools/SAMRecord.h"

/**
 * Represents the results of the reads -> fragment calculation.
 *
 * Contains singleton -- objects whose underlying reads do not overlap their mate pair
 * Contains overlappingPairs -- objects whose underlying reads do overlap their mate pair
 */
template <typename T>
class FragmentCollection {
private:
    std::vector<std::shared_ptr<T>> singletons;
    std::vector<std::pair<std::shared_ptr<T>, std::shared_ptr<T>>> overlappingPairs;

public:
    /**
     * Makes a new collection.
     * Note: this collection stores live pointers to the argument collections.
     * The callers must not modify those arguments after handing them off to this collection.
     *
     * The constructor is private - use the factory method if you need an object.
     */
    FragmentCollection(std::vector<std::shared_ptr<T>>& singletons, std::vector<std::pair<std::shared_ptr<T>, std::shared_ptr<T>>>& overlappingPairs)
    {
        this->singletons = singletons;
        this->overlappingPairs = overlappingPairs;
    }

    /**
     * Create a FragmentCollection containing SAMRecords from a list of reads
     *
     * @param reads a non-null list of reads, ordered by their start location
     * @return a non-null FragmentCollection
     */
    static FragmentCollection<SAMRecord>* create(std::vector<std::shared_ptr<SAMRecord>>& reads)
    {
        int lastStart = -1;
        std::vector<std::shared_ptr<T>> singletons;
        std::vector<std::pair<std::shared_ptr<T>, std::shared_ptr<T>>> overlapping;
        std::unordered_map<std::string, std::shared_ptr<T>> nameMap;

        // build an initial map, grabbing all of the multi-read fragments
        for(auto & read : reads)
        {
            assert(read->getStart() >= lastStart);
            lastStart = read->getStart();

            if(!read->isPaired() || read->mateIsUnmapped() || read->getMateStart() == 0 || read->getMateStart() > read->getEnd())
            {
                // if we know that this read won't overlap its mate, or doesn't have one, jump out early
                singletons.template emplace_back(read);
            } else {
                // the read might overlap it's mate, or is the rightmost read of a pair
                std::string readName = read->getName();
                if(nameMap.template find(readName) != nameMap.end())
                {
                    overlapping.template emplace_back(nameMap.at(readName), read);
                    nameMap.erase(readName);
                } else {
                    nameMap.template emplace(readName, read);
                }
            }
        }

        // add all of the reads that are potentially overlapping but whose mate never showed
        // up to the oneReadPile
        if(!nameMap.empty())
        {
            for(auto iter = nameMap.begin(); iter != nameMap.end(); iter++)
            {
                singletons.template emplace_back(iter->second);
            }
        }

        return new FragmentCollection<SAMRecord>(singletons, overlapping);
    }

    /**
     * Gets the T elements containing overlapping elements, in no particular order
     * The returned collection is unmodifiable.
     */
    std::vector<std::pair<std::shared_ptr<T>, std::shared_ptr<T>>>& getOverlappingPairs(){
        return overlappingPairs;
    }

};


#endif //MUTECT2CPP_MASTER_FRAGMENTCOLLECTION_H
