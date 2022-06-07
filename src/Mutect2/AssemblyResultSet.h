//
// Created by 梦想家xixi on 2021/12/1.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H
#define MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H

#include <map>
#include "AssemblyResult.h"
#include "Haplotype.h"
#include "AssemblyRegion.h"
#include "VariantContext.h"

struct HaplotypeComp {
public:
	bool operator()(const std::shared_ptr<Haplotype> &left, const std::shared_ptr<Haplotype> &right) const {
		if (left->getLength() != right->getLength())
			return left->getLength() > right->getLength();
		int len = left->getLength();
		std::shared_ptr<uint8_t[]> bases1 = left->getBases();
		std::shared_ptr<uint8_t[]> bases2 = right->getBases();
		for (int i = 0; i < len; ++i) {
			if (bases1[i] != bases2[i])
				return bases1[i] < bases2[i];
		}
		return false;
	}
};

struct hash_Haplotype {
	size_t operator()(const std::shared_ptr<Haplotype> &haplotype) const {
		return xxh::xxhash3<64>(haplotype->getBases().get(), haplotype->getBasesLength());
	}
};

struct equal_Haplotype {
	bool operator()(const std::shared_ptr<Haplotype> &left, const std::shared_ptr<Haplotype> &right) const {
		if (left->getLength() != right->getLength())
			return false;
		int size = left->getLength();
		uint8_t *left_bases = left->getBases().get();
		uint8_t *right_bases = right->getBases().get();
		for (int i = 0; i < size; i++) {
			if (left_bases[i] != right_bases[i])
				return false;
		}
		return true;
	}
};

class AssemblyResultSet {
private:
	std::map<int, std::shared_ptr<AssemblyResult>> assemblyResultByKmerSize;
	std::set<std::shared_ptr<Haplotype>, HaplotypeComp> haplotypes;
	std::map<std::shared_ptr<Haplotype>, std::shared_ptr<AssemblyResult>, HaplotypeComp> assemblyResultByHaplotype;
	std::shared_ptr<AssemblyRegion> regionForGenotyping;
	std::shared_ptr<uint8_t[]> fullReferenceWithPadding;
	int fullReferenceWithPaddingLength{};
	std::shared_ptr<SimpleInterval> paddedReferenceLoc;
	bool variationPresent{};
	std::shared_ptr<Haplotype> refHaplotype;
	bool wasTrimmed = false;
	int lastMaxMnpDistanceUsed = -1;
	std::set<int> kmerSizes;
	std::set<std::shared_ptr<VariantContext>, VariantContextComparator> variationEvents;

	bool add(const std::shared_ptr<AssemblyResult> &ar);

	void updateReferenceHaplotype(const std::shared_ptr<Haplotype> &newHaplotype);

	std::vector<std::pair<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>>> calculateOriginalByTrimmedHaplotypes(
			const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion);

	static std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>
	trimDownHaplotypes(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion,
	                   const std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> &haplotypeList);

public:
	AssemblyResultSet() = default;

	bool add(const std::shared_ptr<Haplotype> &h, const std::shared_ptr<AssemblyResult> &ar);

	bool add(const std::shared_ptr<Haplotype> &h);

	void setRegionForGenotyping(std::shared_ptr<AssemblyRegion> regionForGenotyping);

	void setFullReferenceWithPadding(std::shared_ptr<uint8_t[]> fullReferenceWithPadding, int length);

	void setPaddedReferenceLoc(const std::shared_ptr<SimpleInterval> &paddedReferenceLoc);

	std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &getVariationEvents(int maxMnpDistance);

	void regenerateVariationEvents(int distance);

	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> getHaplotypeList();

	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> getSortedHaplotypeList();

	bool isisVariationPresent();

	void deleteEventMap();

	~AssemblyResultSet() = default;

	/**
	 * Trims an assembly result set down based on a new set of trimmed haplotypes.
	 *
	 * @param trimmedAssemblyRegion the trimmed down active region.
	 *
	 * @throws NullPointerException if any argument in {@code null} or
	 *      if there are {@code null} entries in {@code originalByTrimmedHaplotypes} for trimmed haplotype keys.
	 * @throws IllegalArgumentException if there is no reference haplotype amongst the trimmed ones.
	 *
	 * @return never {@code null}, a new trimmed assembly result set.
	 */
	std::shared_ptr<AssemblyResultSet> trimTo(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion);

	/**
	 * Returns the current region for genotyping.
	 *
	 * @return might be {@code null}.
	 */
	std::shared_ptr<AssemblyRegion> getRegionForGenotyping();

	void printSortedHaplotypes();
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H
