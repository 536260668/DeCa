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
	std::string commandLine;
	for (int i = 0; i < argc; ++i) {
		if (i != 0)
			commandLine.append(" ");
		commandLine.append(argv[i]);
	}
	return commandLine;
}

void VCFWriter::writeHeader(const std::string &cmdLine, const std::vector<SAMReadGroupRecord> &readGroup,
                            const std::string &normalSample) {
	mOtherMetaData.insert(std::make_pair("source", "Mutect2"));
	if (!cmdLine.empty()) {
		mOtherMetaData.insert(std::make_pair("Mutect2CommandLine", cmdLine));
	}

	// remove default filter and add a warning
	// TODO:filter not removed
	bcf1_t *rec = bcf_init();
	rec->rid = bcf_hdr_name2id(hdr, "1");
	bcf_remove_filter(hdr, rec, bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS"), 1);

	mOtherMetaData.insert(std::make_pair("filtering_status",
	                                     "Warning: unfiltered Mutect 2 calls.  Please run FilterMutectCalls to remove false positives."));

	// process sample names
	for (const auto &rg: readGroup) {
		std::string sampleName = rg.getAttributesNochange()[SAMReadGroupRecord::READ_GROUP_SAMPLE_TAG];
		if (sampleName == normalSample) {
			mOtherMetaData.insert(std::make_pair("normal_sample", sampleName));
		} else {
			mOtherMetaData.insert(std::make_pair("tumor_sample", sampleName));
		}
		sampleNamesInOrder.push_back(sampleName);
	}
	std::sort(sampleNamesInOrder.begin(), sampleNamesInOrder.end());
	for (const std::string &sampleName: sampleNamesInOrder) {
		char c_str[sampleName.length() + 1];
		strcpy(c_str, sampleName.c_str());
		bcf_hdr_add_sample(hdr, c_str);
	}
	if (bcf_hdr_sync(hdr) < 0)
		throw std::invalid_argument("something wrong when bcf_hdr_sync.");

	// process reference dictionary
	for (auto &sequence: refDict.getSequences()) {
		contigMetaData.push_back(
				"contig=<ID=" + sequence.getSequenceName() + ",length=" + std::to_string(sequence.getSequenceLength()) +
				">");
	}

	// merge and sort all meta date
	//TODO : use bcf_update_info_string()
	for (const auto &infoData: mInfoMetaData) {
		metaDataInSortedOrder.push_back("##" + infoData.second);
	}
	for (const auto &formatData: mFormatMetaData) {
		metaDataInSortedOrder.push_back("##" + formatData.second);
	}
	for (const auto &otherData: mOtherMetaData) {
		metaDataInSortedOrder.push_back("##" + otherData.first + "=" + otherData.second);
	}
	metaDataInSortedOrder.emplace_back("##contig"); // represent all contigs
	std::sort(metaDataInSortedOrder.begin(), metaDataInSortedOrder.end());

	for (const auto &data: metaDataInSortedOrder) {
		if (data != "##contig") {
			appendHeader(data);
			continue;
		}
		for (const auto &contigData: contigMetaData) {
			appendHeader("##" + contigData);
		}
	}


	if (vcf_hdr_write(outFile, hdr) < 0)
		throw std::invalid_argument("something wrong when write vcf header.");
}

void VCFWriter::add(VariantContext vc) {

}

void VCFWriter::close() {
	bcf_hdr_destroy(hdr);
	hts_close(outFile);
}
