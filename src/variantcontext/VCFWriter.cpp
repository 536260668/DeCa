//
// Created by hlf on 6/13/22.
//

#include "VCFWriter.h"
#include <utility>

VCFWriter::VCFWriter(const std::string &outPath, SAMSequenceDictionary sequenceDictionary) : outPath(outPath),
                                                                                             refDict(sequenceDictionary) {
	assert(!outPath.empty());

	if (outPath.length() < 4 || outPath.rfind(".vcf") != (outPath.length() - 4))
		throw std::invalid_argument("something wrong with output file name.");

	if (sequenceDictionary.getSequences().empty())
		throw std::invalid_argument("A reference dictionary is required for creating Tribble indices on the fly");

	char outPath_cstr[outPath.length() + 1];
	strcpy(outPath_cstr, outPath.c_str());
	outFile = hts_open(outPath_cstr, "w");
	if (outFile == nullptr)
		throw std::invalid_argument("something wrong when opening output file.");

	hdr = bcf_hdr_init("w");
	if (hdr == nullptr)
		throw std::invalid_argument("something wrong when initializing vcf header.");
}

void VCFWriter::appendHeader_c(char str[]) {
	if (bcf_hdr_append(hdr, str) < 0)
		throw std::invalid_argument("something wrong when appending vcf header.");
}

void VCFWriter::appendHeader(const std::string &str) {
	char c_str[str.length() + 1];
	strcpy(c_str, str.c_str());
	appendHeader_c(c_str);
}

std::string VCFWriter::getCommandLine(int argc, char *argv[]) {
	std::string commandLine = "##Mutect2CommandLine=";
	for (int i = 0; i < argc; ++i) {
		if (i != 0)
			commandLine.append(" ");
		commandLine.append(argv[i]);
	}
	return commandLine;
}

int VCFWriter::writeHeader(const std::string &cmdLine) {
	appendHeader("##source=Mutect2");
	if (!cmdLine.empty()) {
		appendHeader(cmdLine);
	}
	appendHeader(
			"##filtering_status=Warning: unfiltered Mutect 2 calls.  Please run FilterMutectCalls to remove false positives.");

	bcf_hdr_add_sample(hdr, "/home/cluster/Storage4/hlf/test1223/elwg/ts_tumor_dedup.bam");
	//bcf_hdr_sync(hdr);
	return vcf_hdr_write(outFile, hdr);
}

void VCFWriter::add(VariantContext vc) {

}

void VCFWriter::close() {
	hts_close(outFile);
}
