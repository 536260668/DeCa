//
// Created by 梦想家xixi on 2021/12/1.
//

#include "AssemblyResultSet.h"

#include <utility>
#include "param/ParamUtils.h"

bool AssemblyResultSet::add(const std::shared_ptr<Haplotype> &h, const std::shared_ptr<AssemblyResult> &ar) {
	Mutect2Utils::validateArg(h.get(), "input haplotype cannot be null");
	Mutect2Utils::validateArg(ar.get(), "input assembly-result cannot be null");
	Mutect2Utils::validateArg(h->getGenomeLocation().get(), "the haplotype provided must have a genomic location");

	bool assemblyResultAdditionReturn = add(ar);

	if (haplotypes.find(h) != haplotypes.end()) {
		if (assemblyResultByHaplotype.find(h) == assemblyResultByHaplotype.end()) {
			assemblyResultByHaplotype.insert(std::make_pair(h, ar));
			return true;
		}
		if (assemblyResultByHaplotype.at(h) != ar)
			throw std::invalid_argument("there is already a different assembly result for the input haplotype");
		return assemblyResultAdditionReturn;
	}
	haplotypes.insert(h);
	assemblyResultByHaplotype.insert(std::make_pair(h, ar));
	updateReferenceHaplotype(h);
	if (h->getIsNonReference()) {
		variationPresent = true;
	}
	return true;
}

bool AssemblyResultSet::add(const std::shared_ptr<AssemblyResult> &ar) {
	Mutect2Utils::validateArg(ar.get(), "input assembly-result cannot be null");
	int kmerSize = ar->getKmerSize();
	if (assemblyResultByKmerSize.find(kmerSize) != assemblyResultByKmerSize.end()) {
		if (assemblyResultByKmerSize.at(kmerSize) != ar)
			throw std::invalid_argument("a different assembly result with the same kmerSize was already added");
		return false;
	}
	assemblyResultByKmerSize.insert(std::make_pair(kmerSize, ar));
	kmerSizes.insert(kmerSize);
	return true;
}

void AssemblyResultSet::updateReferenceHaplotype(const std::shared_ptr<Haplotype> &newHaplotype) {
	if (!newHaplotype->getIsReference())
		return;
	if (refHaplotype == nullptr) {
		refHaplotype = newHaplotype;
	} else
		throw std::invalid_argument("the assembly-result-set already have a reference haplotype that is different");
}

void AssemblyResultSet::setRegionForGenotyping(std::shared_ptr<AssemblyRegion> regionForGenotyping) {
	this->regionForGenotyping = regionForGenotyping;
}

void AssemblyResultSet::setFullReferenceWithPadding(std::shared_ptr<uint8_t[]> fullReferenceWithPadding, int length) {
	this->fullReferenceWithPadding = fullReferenceWithPadding;
	this->fullReferenceWithPaddingLength = length;
}

void AssemblyResultSet::setPaddedReferenceLoc(const std::shared_ptr<SimpleInterval> &paddedReferenceLoc) {
	this->paddedReferenceLoc = paddedReferenceLoc;
}

bool AssemblyResultSet::add(const std::shared_ptr<Haplotype> &h) {
	Mutect2Utils::validateArg(h.get(), "input haplotype can not be null");
	Mutect2Utils::validateArg(h->getGenomeLocation().get(), "haplotype genomeLocation cannot be null");
	if (haplotypes.find(h) != haplotypes.end())
		return false;
	haplotypes.insert(h);
	updateReferenceHaplotype(h);
	return true;
}

std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &
AssemblyResultSet::getVariationEvents(int maxMnpDistance) {
	ParamUtils::isPositiveOrZero(maxMnpDistance, "maxMnpDistance may not be negative.");
	bool sameMnpDistance = lastMaxMnpDistanceUsed != -1 && maxMnpDistance == lastMaxMnpDistanceUsed;
	lastMaxMnpDistanceUsed = maxMnpDistance;
	bool flag = false;
	for (const auto &haplotype: haplotypes) {
		if (haplotype->getIsNonReference() &&
		    (haplotype->getEventMap() == nullptr || haplotype->getEventMap()->empty())) {
			flag = true;
			break;
		}
	}
	if (variationEvents.empty() || !sameMnpDistance || flag) {
		regenerateVariationEvents(maxMnpDistance);
	}
	return variationEvents;
}

void AssemblyResultSet::regenerateVariationEvents(int maxMnpDistance) {
	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> haplotypeList = getHaplotypeList();
	EventMap::buildEventMapsForHaplotypes(*haplotypeList, fullReferenceWithPadding, fullReferenceWithPaddingLength,
	                                      paddedReferenceLoc, false, maxMnpDistance);
	variationEvents = EventMap::getAllVariantContexts(*haplotypeList);
	lastMaxMnpDistanceUsed = maxMnpDistance;
	for (const std::shared_ptr<Haplotype> &haplotype: *haplotypeList) {
		if (haplotype->getIsNonReference()) {
			variationPresent = true;
			break;
		}
	}
}

std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>>
AssemblyResultSet::getHaplotypeList() {
	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>>
			res = std::make_shared<std::vector<std::shared_ptr<Haplotype>>>();
	res->reserve(haplotypes.size());
	for (const auto &haplotype: haplotypes) {
		res->emplace_back(haplotype);
	}
	return res;
}

std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>
AssemblyResultSet::calculateOriginalByTrimmedHaplotypes(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion) {
	const std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> haplotypeList = getHaplotypeList();
	std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype >> originalByTrimmedHaplotypes = trimDownHaplotypes(
			trimmedAssemblyRegion, haplotypeList);
	std::vector<std::shared_ptr<Haplotype>> trimmedHaplotypes;
	trimmedHaplotypes.reserve(originalByTrimmedHaplotypes->size());
	for (const auto &element: *originalByTrimmedHaplotypes) {
		trimmedHaplotypes.emplace_back(element.first);
	}
	std::sort(trimmedHaplotypes.begin(), trimmedHaplotypes.end(),
	          [](const std::shared_ptr<Haplotype> &left, const std::shared_ptr<Haplotype> &right) {
		          return *left < *right;
	          });
	std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype >>
			sortedOriginalByTrimmedHaplotypes = mapOriginalToTrimmed(originalByTrimmedHaplotypes, trimmedHaplotypes);
	return sortedOriginalByTrimmedHaplotypes;
}

std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>
AssemblyResultSet::trimDownHaplotypes(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion,
                                      const std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> &haplotypeList) {
	std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>
			originalByTrimmedHaplotypes = std::make_shared<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>();
	//todo: memory leak here?
	for (const auto &h: *haplotypeList) {
		std::shared_ptr<Haplotype> trimmed = h->trim(trimmedAssemblyRegion->getExtendedSpan());
		if (trimmed != nullptr) {
			auto iter = originalByTrimmedHaplotypes->find(trimmed);
			if (iter != originalByTrimmedHaplotypes->end()) {
				if (trimmed->getIsReference()) {
					originalByTrimmedHaplotypes->erase(iter);
					originalByTrimmedHaplotypes->insert({trimmed, h});
				}
			} else {
			    originalByTrimmedHaplotypes->insert({trimmed, h});
			}
		} else if (h->getIsReference())
			throw std::invalid_argument("trimming eliminates the reference haplotype");
	}
	return originalByTrimmedHaplotypes;
}

std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>
AssemblyResultSet::mapOriginalToTrimmed(
		const std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>> &originalByTrimmedHaplotypes,
		const std::vector<std::shared_ptr<Haplotype>> &trimmedHaplotypes) {
	std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>
			sortedOriginalByTrimmedHaplotypes = std::make_shared<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>();
	sortedOriginalByTrimmedHaplotypes->reserve(trimmedHaplotypes.size());
	for (const auto &trimmed: trimmedHaplotypes) {
		sortedOriginalByTrimmedHaplotypes->insert({trimmed, originalByTrimmedHaplotypes->at(trimmed)});
	}
	return sortedOriginalByTrimmedHaplotypes;
}

std::shared_ptr<AssemblyResultSet>
AssemblyResultSet::trimTo(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion) {
	std::shared_ptr<std::unordered_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype >>
			originalByTrimmedHaplotypes = calculateOriginalByTrimmedHaplotypes(trimmedAssemblyRegion);
	if (refHaplotype == nullptr)
		throw std::invalid_argument("refHaplotype is null");

	std::shared_ptr<AssemblyResultSet> result = std::make_shared<AssemblyResultSet>();
	for (const auto &element: *originalByTrimmedHaplotypes) {
		const std::shared_ptr<Haplotype> &trimmed = element.first;
		const std::shared_ptr<Haplotype> &original = element.second;
		if (original == nullptr)
			throw std::invalid_argument("all trimmed haplotypes must have an original one");

		if(assemblyResultByHaplotype.find(original) == assemblyResultByHaplotype.end())
		    result->add(trimmed);
		else {
		    std::shared_ptr<AssemblyResult> &as = assemblyResultByHaplotype.at(original);
		    result->add(trimmed, as);
		}
	}

	result->setRegionForGenotyping(trimmedAssemblyRegion);
	result->setFullReferenceWithPadding(fullReferenceWithPadding, fullReferenceWithPaddingLength);
	result->setPaddedReferenceLoc(paddedReferenceLoc);
	result->variationPresent = false;
	for (const auto &haplotype: haplotypes) {
		if (haplotype->getIsNonReference()) {
			result->variationPresent = true;
		}
	}
	result->wasTrimmed = true;
	return result;
}

bool AssemblyResultSet::isisVariationPresent() {
	return variationPresent && haplotypes.size() > 1;
}

std::shared_ptr<AssemblyRegion> AssemblyResultSet::getRegionForGenotyping() {
	return regionForGenotyping;
}

void AssemblyResultSet::deleteEventMap() {
	auto haplotypesToReleased = *this->getHaplotypeList();
	for (auto &item: haplotypesToReleased) {
		if (item->getEventMap() != nullptr){
			delete item->getEventMap();
			item->setEventMap(nullptr);
		}
	}
}
