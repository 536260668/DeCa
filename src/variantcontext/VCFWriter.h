//
// Created by hlf on 6/13/22.
//

#ifndef MUTECT2CPP_MASTER_VCFWRITER_H
#define MUTECT2CPP_MASTER_VCFWRITER_H

#include "htslib/vcf.h"
#include "VariantContext.h"
#include <string>

class VCFWriter {
	std::string VERSION_LINE;

	public:
	void add(VariantContext vc);
};

#endif //MUTECT2CPP_MASTER_VCFWRITER_H
