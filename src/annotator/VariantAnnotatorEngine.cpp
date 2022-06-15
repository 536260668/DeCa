//
// Created by lhh on 5/20/22.
//

#include "variantcontext/builder/VariantContextBuilder.h"
#include "VariantAnnotatorEngine.h"

VariantAnnotatorEngine::VariantAnnotatorEngine() {}

VariantAnnotatorEngine::VariantAnnotatorEngine(std::vector<shared_ptr<InfoFieldAnnotation>>&& InfoFieldAnnotationList, std::vector<shared_ptr<GenotypeAnnotation>>&& GenotypeAnnotationList): infoAnnotations(InfoFieldAnnotationList), genotypeAnnotations(GenotypeAnnotationList),
useRawAnnotations(false), keepRawCombinedAnnotations(false)
{

}

shared_ptr<VariantContext> VariantAnnotatorEngine::annotateContext(shared_ptr<VariantContext> vc, ReferenceContext& ref,
                                                                    AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    assert(vc != nullptr);

    // annotate genotypes, creating another new VC in the process
    VariantContextBuilder builder(vc);
    builder.setGenotypes(annotateGenotypes(ref, vc, likelihoods));
    auto newGenotypeAnnotatedVC = builder.make();

    auto infoAnnotMap = make_shared<std::map<std::string, AttributeValue>>();

    for(auto& annotationType : infoAnnotations)
    {
        auto annotationsFromCurrentType = annotationType->annotate(make_shared<ReferenceContext>(ref), newGenotypeAnnotatedVC, likelihoods);
        if(annotationsFromCurrentType != nullptr)
        {
            for(auto& annotation : *annotationsFromCurrentType)
            {
                infoAnnotMap->emplace(annotation);
            }
        }
    }
    auto annotated = builder.setAttributes(infoAnnotMap)->make();

    return annotated;
}

GenoTypesContext *VariantAnnotatorEngine::annotateGenotypes(ReferenceContext &ref, shared_ptr<VariantContext> vc,
                                                             AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    if(genotypeAnnotations.empty())
        return vc->getGenotypes();

    GenoTypesContext* genotypes = new GenoTypesContext(vc->getNSamples());
    for(int i=0; i<genotypes->getSize(); i++)
    {
        auto genotype = genotypes->get(i);
        GenotypeBuilder gb(genotype);
        for(auto& annotation  : genotypeAnnotations)
        {
            annotation->annotate(ref, vc, genotype, gb, likelihoods);
        }
        genotypes->add(gb.make());
    }

    return genotypes;
}