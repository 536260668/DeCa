//
// Created by lhh on 10/25/21.
//

#include "BandPassActivityProfile.h"
#include "MathUtils.h"

BandPassActivityProfile::BandPassActivityProfile(int maxProbPropagationDistance, double activeProbThreshold,
                                                 int maxFilterSize, double sigma, bool adaptiveFilterSize,
                                                 SAMFileHeader *header) : ActivityProfile(maxProbPropagationDistance, activeProbThreshold, header), sigma(sigma)
{
    std::vector<double> * fullKernel = makeKernel(maxFilterSize, sigma);
    filterSize = adaptiveFilterSize ? determineFilterSize(fullKernel, MIN_PROB_TO_KEEP_IN_FILTER) : maxFilterSize;
    gaussianKernel = makeKernel(filterSize, sigma);
    delete fullKernel;
}

std::vector<double> * BandPassActivityProfile::makeKernel(int filterSize, double sigma)
{
    int bandSize = (filterSize << 1) + 1;
    auto * kernel = new std::vector<double>;
    kernel->reserve(bandSize);
    for(int i=0; i<bandSize; i++)
    {
        kernel->push_back(MathUtils::normalDistribution(filterSize, sigma, i));   //---to be tested
    }
    return MathUtils::normalizeSumToZero(kernel);
}

BandPassActivityProfile::~BandPassActivityProfile()
{
    delete gaussianKernel;  //---deleting null pointer has no effect
}

int BandPassActivityProfile::determineFilterSize(std::vector<double> * kernel, double minProbToKeepInFilter)
{
    int middle = (kernel->size() - 1) >> 1;
    int filterEnd = middle;
    while (filterEnd > 0){
        if(kernel->operator[](filterEnd - 1) < minProbToKeepInFilter)
            break;
        filterEnd--;
    }
    return middle - filterEnd;
}

unique_ptr<vector<struct ActivityProfileState>> BandPassActivityProfile::processState(ActivityProfileState &justAddedState)
{
    unique_ptr<vector<ActivityProfileState>> states(new std::vector<ActivityProfileState>);

    double activeProb = justAddedState.isActiveProb();

    if(activeProb > 0.0)
    {
        for(int i = -filterSize; i <= filterSize; i++)
        {
            optional<SimpleInterval> loc = getLocForOffset(justAddedState.getLoc(), i);
            if(loc.has_value())
            {
                double newProb = activeProb * gaussianKernel->operator[](i + filterSize);
                states->emplace_back(ActivityProfileState(loc.value(), newProb));
            }
        }
    } else{
        states->emplace_back(justAddedState);
    }
    return states;
}

int BandPassActivityProfile::getMaxProbPropagationDistance()
{
    return maxProbPropagationDistance + filterSize;
}