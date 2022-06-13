//
// Created by lhh on 6/2/22.
//

#ifndef MUTECT2CPP_MASTER_BASEQUALITY_H
#define MUTECT2CPP_MASTER_BASEQUALITY_H

#include "PreAlleleAnnotation.h"

/**
 * Median base quality of bases supporting each allele.
 *
 * <p>The output is an array containing, for each allele, the median base quality at the variant (for SNVs) and one base before the variant (for indels) over all reads that best match that allele.</p>
 *
 * <p>For example, for variant context with ref allele A and alt allele C we use base qualities at the A and C.  For variant context with ref allele AG and alt allele A (deletion),
 * we use base qualities at the A.  For variant context with ref allele A and alt allele AG (insertion) we use base qualities at the A.</p>
 */
class BaseQuality : public PreAlleleAnnotation{
public:

};


#endif //MUTECT2CPP_MASTER_BASEQUALITY_H
