//
// Created by 梦想家xixi on 2021/10/30.
//

#include "QualityUtils.h"
#include "Mutect2Utils.h"
#include <cmath>

double QualityUtils::qualToErrorProbCache[255] {0};

uint8_t QualityUtils::errorProbToQual(double prob, uint8_t maxQual) {
    Mutect2Utils::validateArg(Mutect2Utils::goodProbability(prob), "errorRate must be good probability");
    double d = std::round(-10.0*std::log10(prob));
    return boundQual((int)d, maxQual);
}

uint8_t QualityUtils::boundQual(int qual, uint8_t maxQual) {
    return (uint8_t)(std::max(std::min(qual, maxQual & 0xff), 1) & 0xff);
}

uint8_t QualityUtils::errorProbToQual(double errorRate) {
    return errorProbToQual(errorRate, 93);
}

void QualityUtils::initial() {
    for(int i = 0; i < 255; i++) {
        qualToErrorProbCache[i] = qualToErrorProb((double) i);
    }
}

double QualityUtils::qualToErrorProb(double qual) {
    Mutect2Utils::validateArg(qual >= 0.0, "Qual must be >= 0.0");
    return std::pow(10.0, qual / -10.0);
}

double QualityUtils::qualToErrorProb(uint8_t qual) {
    return qualToErrorProbCache[(int)qual & 0xff];
}
