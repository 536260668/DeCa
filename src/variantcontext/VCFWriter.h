//
// Created by hlf on 6/13/22.
//

#ifndef MUTECT2CPP_MASTER_VCFWRITER_H
#define MUTECT2CPP_MASTER_VCFWRITER_H

#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "VariantContext.h"
#include "samtools/SAMSequenceDictionary.h"
#include "VCFConstants.h"
#include <string>

class VCFWriter {
	std::string VERSION_LINE;
	std::string outPath;
	SAMSequenceDictionary refDict;
	htsFile *outFile;
	bcf_hdr_t *hdr;

public:

	VCFWriter(const std::string &outPath, SAMSequenceDictionary sequenceDictionary);

	static std::string getCommandLine(int argc, char *argv[]);

	int writeHeader(const std::string &cmdLine);

	void add(VariantContext vc);

	void close();

private:

	void appendHeader_c(char str[]);

	void appendHeader(const std::string &str);
};

#endif //MUTECT2CPP_MASTER_VCFWRITER_H
