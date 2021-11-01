//
// Created by lhh on 10/25/21.
//

#include <limits>
#include "ActivityProfile.h"
#include "assert.h"


ActivityProfile::ActivityProfile(int maxProbPropagationDistance, double activeProbThreshold, sam_hdr_t * header):
                    stateList(), maxProbPropagationDistance(maxProbPropagationDistance),
                    activeProbThreshold(activeProbThreshold) , regionStartLoc(), regionStopLoc() ,samHeader(header),
                    contigLength(-1), regions(new vector<AssemblyRegion>)
{
}

ActivityProfile::~ActivityProfile()
{
    delete regions;
}

bool ActivityProfile::isEmpty()
{
    return stateList.empty();
}

void ActivityProfile::add(ActivityProfileState & state)
{

    if(regionStartLoc.getContig().empty()){
        regionStartLoc = state.getLoc();
        regionStopLoc = regionStartLoc;
        contigLength = sam_hdr_tid2len(samHeader, sam_hdr_name2tid(samHeader, regionStartLoc.getContig().c_str()));
    } else {
        // GATK requires the activityProfile to be contiguous with the previously added one
        assert(regionStopLoc.getStart() == state.getLoc().getStart() - 1);
        regionStopLoc = state.getLoc();
    }

    unique_ptr<vector<ActivityProfileState>> processedStates = processState(state);
    for(ActivityProfileState & processedState : *processedStates){
        incorporateSingleState(processedState);
    }
}

unique_ptr<vector<ActivityProfileState>> ActivityProfile::processState(ActivityProfileState & justAddedState)
{

}

optional<SimpleInterval> ActivityProfile::getLocForOffset(SimpleInterval relativeLoc, int offset)
{
    int start = relativeLoc.getStart() + offset;
    if(start < 1 || start > contigLength)
        return std::nullopt;
    else
        return std::optional<SimpleInterval>{SimpleInterval(regionStartLoc.getContig(), start, start)};

}

void ActivityProfile::incorporateSingleState(ActivityProfileState &stateToAdd)
{
    int size = stateList.size();
    int position = stateToAdd.getOffset(&regionStartLoc);
    //assert(position <= size);

    if(position >= 0){
        if(position < size)
        {
            stateList[position].setIsActiveProb(stateList[position].isActiveProb() + stateToAdd.isActiveProb());
        } else{
            stateList.emplace_back(stateToAdd);
        }

    }
}

int ActivityProfile::getEnd()
{
    return regionStopLoc.getEnd();
}

vector<AssemblyRegion> * ActivityProfile::popReadyAssemblyRegions(int assemblyRegionExtension, int minRegionSize,
                                                                  int maxRegionSize, bool forceConversion)
{
    if(!regions->empty())
    {
        regions->clear();
    }

    while(true)
    {
        optional<AssemblyRegion> ReadyAssemblyRegion = popReadyAssemblyRegion(assemblyRegionExtension, minRegionSize, maxRegionSize, forceConversion);
        if(ReadyAssemblyRegion.has_value()){
            // only active regions are added here
            if(ReadyAssemblyRegion.value().getIsActive())
                regions->emplace_back(ReadyAssemblyRegion.value());
            cout << ReadyAssemblyRegion.value();
        } else {
            return regions;
        }
    }
}

optional<struct AssemblyRegion> ActivityProfile::popReadyAssemblyRegion(int assemblyRegionExtension, int minRegionSize, int maxRegionSize,
                                        bool forceConversion)
{
    if(stateList.empty()){
        return nullopt;
    }

    ActivityProfileState first = stateList[0];
    bool isActiveRegion = first.isActiveProb() > activeProbThreshold;
    int offsetOfNextRegionEnd = findEndOfRegion(isActiveRegion, minRegionSize, maxRegionSize, forceConversion);
    if(offsetOfNextRegionEnd == -1){
        return nullopt;
    }

    vector<ActivityProfileState> sub = {stateList.begin(), stateList.begin() + offsetOfNextRegionEnd + 1};
    stateList = {stateList.begin() + offsetOfNextRegionEnd + 1, stateList.end()};

    if(stateList.empty())
    {
        regionStartLoc.clearContig();
        regionStopLoc.clearContig();
        assert(regionStartLoc.getContig().empty());

    } else {
        regionStartLoc = stateList[0].getLoc();
    }


    SimpleInterval regionLoc = SimpleInterval(first.getLoc().getContig(), first.getLoc().getStart(), first.getLoc().getStart() + offsetOfNextRegionEnd);
    return optional<AssemblyRegion>( AssemblyRegion(regionLoc, sub, isActiveRegion, assemblyRegionExtension, samHeader));
}

int ActivityProfile::findEndOfRegion(bool isActiveRegion, int minRegionSize, int maxRegionSize, bool forceConversion)
{
    if(!forceConversion && stateList.size() < maxRegionSize + getMaxProbPropagationDistance())
        return -1;

    int endOfActiveRegion = findFirstActivityBoundary(isActiveRegion, maxRegionSize);

    if(isActiveRegion && endOfActiveRegion == maxRegionSize){
        endOfActiveRegion = findBestCutSite(endOfActiveRegion, minRegionSize);
    }

    return endOfActiveRegion - 1;
}

int ActivityProfile::findFirstActivityBoundary(bool isActiveRegion, int maxRegionSize)
{
    int nStates = stateList.size();
    int endOfActiveRegion = 0;

    while (endOfActiveRegion < nStates && endOfActiveRegion < maxRegionSize){
        if(stateList[endOfActiveRegion].isActiveProb() > activeProbThreshold != isActiveRegion)
            break;
        endOfActiveRegion++;
    }

    return endOfActiveRegion;
}

int ActivityProfile::findBestCutSite(int endOfActiveRegion, int minRegionSize)
{
    assert(endOfActiveRegion >= minRegionSize);
    assert(minRegionSize >= 0);

    int minI = endOfActiveRegion - 1;
    double minP = numeric_limits<double>::max();
    for(int i = minI; i >= minRegionSize - 1; i--)
    {
        double cur = stateList[i].isActiveProb();
        if(cur < minP && isMinimum(i)){
            minP = cur;
            minI = i;
        }
    }

    return minI + 1;
}

bool ActivityProfile::isMinimum(int index)
{
    if(index == stateList.size() - 1)
        // we cannot be at a minimum if the current position is the last in the state list
        return false;
    else if(index < 1)
        // we cannot be at a minimum if the current position is the first or second
        return false;
    else{
        double indexP = stateList[index].isActiveProb();
        return indexP <= stateList[index+1].isActiveProb() && indexP < stateList[index-1].isActiveProb();
    }
}