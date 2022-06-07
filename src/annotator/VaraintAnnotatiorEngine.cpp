//
// Created by lhh on 5/20/22.
//

#include "variantcontext/builder/VariantContextBuilder.h"
#include "VaraintAnnotatiorEngine.h"

VaraintAnnotatiorEngine::VaraintAnnotatiorEngine() {}

shared_ptr<VariantContext> VaraintAnnotatiorEngine::annotateContext(shared_ptr<VariantContext> vc, ReferenceContext& ref,
                                                                    AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    assert(vc != nullptr);

    // annotate genotypes, creating another new VC in the process
    VariantContextBuilder builder(vc);
    // TODO: finish this method 2022.6.2

    return vc;
}

GenoTypesContext *VaraintAnnotatiorEngine::annotateGenotypes(ReferenceContext &ref, shared_ptr<VariantContext> vc,
                                                             AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    if(genotypeAnnotations.empty())
        return vc->getGenotypes();

    GenoTypesContext* genotypes = new GenoTypesContext(vc->getNSamples());
    for(int i=0; i<genotypes->getSize(); i++)
    {
        auto genotype = genotypes->get(i);
        GenotypeBuilder gb(genotype);
        for(GenotypeAnnotation* annotation : genotypeAnnotations)
        {
            annotation->annotate(ref, vc, genotype, gb, likelihoods);
        }
        genotypes->add(gb.make());
    }

    return genotypes;
}