//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_BASEQUALITYFILTER_H
#define MUTECT2CPP_MASTER_BASEQUALITYFILTER_H

#include "HardFilter.h"


class BaseQualityFilter : public HardFilter{
private:
    double minMedianBaseQuality;

public:
    BaseQualityFilter(double minMedianBaseQuality);
    virtual ErrorType errorType();
    virtual bool isArtifact(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine);
    virtual std::vector<std::string> requiredAnnotations();
    virtual std::string filterName();
};


#endif //MUTECT2CPP_MASTER_BASEQUALITYFILTER_H
