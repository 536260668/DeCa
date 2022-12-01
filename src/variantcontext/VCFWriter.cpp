//
// Created by hlf on 6/13/22.
//

#include "VCFWriter.h"
#include <utility>
#include "VCFConstants.h"

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


void VCFWriter::add(std::shared_ptr<VariantContext>& vc) {
    bcf1_t * hts_vc = bcf_init();
    hts_vc->pos = vc->getStart();
    hts_vc->rid = vc->getContigInt();
    hts_vc->n_sample = bcf_hdr_nsamples(hdr);
    bcf_float_set_missing(hts_vc->qual);
    hts_vc->rlen = vc->getReference()->getLength();
    int len = 0;
    len += (vc->getReference()->getLength()+1);
    for(auto allele : vc->getAlternateAlleles()) {
        len += (allele->getLength()+1);
    }
    char * tmp = new char [len]{0};
    int sum = 0;
    memcpy(tmp, vc->getReference()->getBases().get(), vc->getReference()->getLength());
    sum += vc->getReference()->getLength();
    tmp[sum] = ',';
    sum++;
    for(auto allele : vc->getAlternateAlleles()) {
        memcpy(tmp + sum, allele->getBases().get(), allele->getLength());
        sum += allele->getLength();
        tmp[sum] = ',';
        sum++;
    }
    tmp[sum-1] = '\0';
    bcf_update_alleles_str(hdr, hts_vc, tmp);

    //
    for(auto & kv : vc->getAttributes()) {
        if(kv.first == VCFConstants::DEPTH_KEY) {
            int x = kv.second.getAttributeAsInt();
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::DEPTH_KEY.c_str(), &x, 1);
            continue;
        }
        if(kv.first == VCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY) {
            int x = kv.second.getAttributeAsInt();
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY.c_str(), &x, 1);
            continue;
        }
        if(kv.first == VCFConstants::MEDIAN_READ_POSITON_KEY) {
            int x = kv.second.getAttributeAsInt();
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::MEDIAN_READ_POSITON_KEY.c_str(), &x, 1);
            continue;
        }
        if(kv.first == VCFConstants::MEDIAN_BASE_QUALITY_KEY) {
            std::vector<int> x = kv.second.getAttributeAsIntVector();
            int list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::MEDIAN_BASE_QUALITY_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY) {
            std::vector<int> x = kv.second.getAttributeAsIntVector();
            int list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::MEDIAN_MAPPING_QUALITY_KEY) {
            std::vector<int> x = kv.second.getAttributeAsIntVector();
            int list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::MEDIAN_MAPPING_QUALITY_KEY.c_str(), list, x.size());
            continue;
        }
    }
    bcf_add_filter(hdr, hts_vc, 0);

    bcf_write(outFile, hdr, hts_vc);
    bcf_destroy(hts_vc);
    delete[] tmp;
}

void VCFWriter::close() {
	bcf_hdr_destroy(hdr);
	hts_close(outFile);
}

void
VCFWriter::writeHeader(const string &cmdLine, const vector<SAMReadGroupRecord> &readGroup, const string &normalSample) {
    //version
    std::string version = "VCFv4.2";
    bcf_hdr_set_version(hdr, version.c_str());

    //format
    for(auto & val : mFormatMetaData) {
        string format = "##" + val.second;
        int len = format.size();
        bcf_hrec_t * bht = bcf_hdr_parse_line(hdr, format.c_str(), &len);
        bcf_hdr_add_hrec(hdr, bht);
    }

    //info
    for(auto & val : mInfoMetaData) {
        string info = "##" + val.second;
        int len = info.size();
        bcf_hrec_t * bht = bcf_hdr_parse_line(hdr, info.c_str(), &len);
        bcf_hdr_add_hrec(hdr, bht);
    }

    //contig
    for (auto &sequence: refDict.getSequences()) {
        contigMetaData.push_back("contig=<ID=" + sequence.getSequenceName() + ",length=" + std::to_string(sequence.getSequenceLength()) +">");
	}
    for(auto & val : contigMetaData) {
        string contig = "##" + val;
        int len = contig.size();
        bcf_hrec_t * bht = bcf_hdr_parse_line(hdr, contig.c_str(), &len);
        bcf_hdr_add_hrec(hdr, bht);
    }

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
    for(auto & k : mOtherMetaData) {
        std::string sample = "##" + k.first + "=" + k.second;
        int len = sample.size();
        bcf_hrec_t * bht = bcf_hdr_parse_line(hdr, sample.c_str(), &len);
        bcf_hdr_add_hrec(hdr, bht);
    }

    bcf_hdr_sync(hdr);
    vcf_hdr_write(outFile, hdr);
}
