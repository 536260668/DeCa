//
// Created by lhh on 6/15/22.
//

#include "TandemRepeat.h"
#include "VCFConstants.h"
#include "utils/variant/GATKVariantContextUtils.h"

std::vector<std::string> TandemRepeat::getKeyNames() {
    return {VCFConstants::STR_PRESENT_KEY,
            VCFConstants::REPEAT_UNIT_KEY,
            VCFConstants::REPEATS_PER_ALLELE_KEY};
}

std::shared_ptr<std::map<std::string, AttributeValue>>
TandemRepeat::annotate(shared_ptr<ReferenceContext> ref, shared_ptr<VariantContext> vc,
                       AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    //todo::filter
    assert(vc != nullptr);
    if ( !vc->isIndel()) {
        return nullptr;
    }

    return nullptr;
}

