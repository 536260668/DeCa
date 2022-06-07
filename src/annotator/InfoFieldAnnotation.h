//
// Created by lhh on 5/20/22.
//

#ifndef MUTECT2CPP_MASTER_INFOFIELDANNOTATION_H
#define MUTECT2CPP_MASTER_INFOFIELDANNOTATION_H

#include "VariantAnnotation.h"
#include "engine/ReferenceContext.h"
#include "VariantContext.h"
#include "Genotype.h"
#include "variantcontext/builder/GenotypeBuilder.h"
#include "utils/genotyper/AlleleLikelihoods.h"

/**
 * Annotations relevant to the INFO field of the variant file (ie annotations for sites).
 */
class InfoFieldAnnotation : public VariantAnnotation{
public:
    virtual void annotate(shared_ptr<ReferenceContext> ref, shared_ptr<VariantContext> vc, AlleleLikelihoods<SAMRecord, Allele>* likelihoods) = 0;

    std::string toString();
};


#endif //MUTECT2CPP_MASTER_INFOFIELDANNOTATION_H
