//
// Created by lhh on 6/15/22.
//

#ifndef MUTECT2CPP_MASTER_TANDEMREPEAT_H
#define MUTECT2CPP_MASTER_TANDEMREPEAT_H

#include "InfoFieldAnnotation.h"

class TandemRepeat : public InfoFieldAnnotation{
public:
    std::vector<std::string> getKeyNames();

    std::shared_ptr<std::map<std::string, AttributeValue>> annotate(shared_ptr<ReferenceContext> ref, shared_ptr<VariantContext> vc, AlleleLikelihoods<SAMRecord, Allele>* likelihoods);
};


#endif //MUTECT2CPP_MASTER_TANDEMREPEAT_H
