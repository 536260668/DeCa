//
// Created by hlf on 8/16/22.
//

#ifndef MUTECT2CPP_MASTER_READFORPAIRHMM_H
#define MUTECT2CPP_MASTER_READFORPAIRHMM_H

#include "xxhash.hpp"
#include "parallel_hashmap/phmap.h"

struct ReadForPairHMM {
	int rslen;
	const uint8_t *q, *i, *d, *c;
	const uint8_t *rs;
	shared_ptr<uint32_t[]> charCombination;
	xxh::hash64_t hashCode;

	ReadForPairHMM(int _rslen, const uint8_t *readQuals, const uint8_t *insGops, const uint8_t *delGops,
	               const char *gapConts, const uint8_t *reads) : rslen(_rslen), q(readQuals), i(insGops), d(delGops),
	                                                             c((uint8_t *) gapConts), rs(reads) {
		/*
		 * Calculate the hashcode and store it
		 * The calculation of hashcode does not consider the base sequence,
		 * which can remove most of the duplication and improve the efficiency.
		 * But the base sequence should be considered when determining whether the two are equal.
		 * */

		// 4 * sizeof(uint8_t) * rslen = sizeof(uint32_t) * rslen
		charCombination = std::shared_ptr<uint32_t[]>(new uint32_t[rslen]);
		auto *charP = (char *) charCombination.get();
		memcpy(charP, d, rslen);
		memcpy(charP + rslen, i, rslen);
		memcpy(charP + 2 * rslen, c, rslen);
		memcpy(charP + 3 * rslen, q, rslen);
		for (int k = 0; k < rslen; ++k) {
			charCombination[k] &= 2139062143;   // means every char &= 127
		}
		hashCode = xxh::xxhash3<64>(charCombination.get(), rslen);
	}
};

struct ReadForPairHMMHash {
	xxh::hash64_t operator()(const std::shared_ptr<ReadForPairHMM> &t) const {
		return t->hashCode;
	}
};

struct ReadForPairHMMEqual {
	bool operator()(const std::shared_ptr<ReadForPairHMM> &t1, const std::shared_ptr<ReadForPairHMM> &t2) const {
		if (t1->rslen != t2->rslen)
			return false;

		int len = t1->rslen;

		if (memcmp(t1->rs, t2->rs, len) != 0)
			return false;

		if (memcmp(t1->charCombination.get(), t2->charCombination.get(), 4 * len) != 0)
			for (int i = 0; i < len; ++i)
				if ((t1->charCombination[i] & 2139062143) != (t2->charCombination[i] & 2139062143))
					return false;

		return true;
	}
};


#endif //MUTECT2CPP_MASTER_READFORPAIRHMM_H
