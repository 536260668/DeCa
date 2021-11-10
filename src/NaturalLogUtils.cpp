//
// Created by 梦想家xixi on 2021/11/1.
//

#include "NaturalLogUtils.h"
#include "Mutect2Utils.h"

double NaturalLogUtils::qualToLogProbCache[255] {0};

double NaturalLogUtils::qualToLogErrorProb(double qual) {
    Mutect2Utils::validateArg(qual >= 0.0, "qual must be >= 0.0 ");
    return qual * PHRED_TO_LOG_ERROR_PROB_FACTOR;
}

double NaturalLogUtils::qualToLogErrorProb(uint8_t qual) {
    return qualToLogErrorProb((double)(qual & 0xff));
}

double NaturalLogUtils::log1mexp(double a) {
    if (a > 0) return std::nan("");
    if (a == 0) return -std::numeric_limits<double>::infinity();
    return (a < LOG1MEXP_THRESHOLD) ? std::log1p(-std::exp(a)) : std::log(-std::expm1(a));
}

void NaturalLogUtils::initial() {
    for(int i = 0; i < 255; i++) {
        qualToLogProbCache[i] = log1mexp(qualToLogErrorProb(double(i)));
    }
}

double NaturalLogUtils::qualToLogProb(uint8_t qual) {
    return qualToLogProbCache[(int) qual & 0xff];
}
