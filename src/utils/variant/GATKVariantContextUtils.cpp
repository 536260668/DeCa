//
// Created by lhh on 4/26/22.
//

#include <cassert>
#include <utils/Utils.h>
#include "GATKVariantContextUtils.h"

int GATKVariantContextUtils::findNumberOfRepetitions(uint8_t *repeatUnitFull, int repeatUnitFullLength, int offsetInRepeatUnitFull,
                                                     int repeatUnitLength, uint8_t *testStringFull, int testStringFullLength,
                                                     int offsetInTestStringFull, int testStringLength,
                                                     bool leadingRepeats) {
    if (testStringLength == 0){
        return 0;
    }

    assert(repeatUnitLength >= 0 && repeatUnitLength <= repeatUnitFullLength);
    assert(offsetInRepeatUnitFull >= 0 && offsetInRepeatUnitFull < repeatUnitFullLength);
    assert(offsetInTestStringFull >= 0 && offsetInTestStringFull < testStringFullLength);
    assert(testStringLength >= 0 && testStringLength <= testStringFullLength);

    int lengthDifference = testStringLength - repeatUnitLength;

    if(leadingRepeats)
    {
        int numRepeats = 0;
        // look forward on the test string
        for (int start = 0; start <= lengthDifference; start += repeatUnitLength) {
            if(Utils::equalRange(testStringFull, start + offsetInTestStringFull, repeatUnitFull, offsetInRepeatUnitFull, repeatUnitLength)) {
                numRepeats++;
            } else {
                return numRepeats;
            }
        }
        return numRepeats;
    } else {
        // look backward. For example, if repeatUnit = AT and testString = GATAT, number of repeat units is still 2
        int numRepeats = 0;
        // look backward on the test string
        for (int start = lengthDifference; start >= 0; start -= repeatUnitLength) {
            if (Utils::equalRange(testStringFull, start + offsetInTestStringFull, repeatUnitFull, offsetInRepeatUnitFull, repeatUnitLength)) {
                numRepeats++;
            } else {
                return numRepeats;
            }
        }
        return numRepeats;
    }

}

int GATKVariantContextUtils::findNumberOfRepetitions(uint8_t* repeatUnit, int repeatUnitLength, uint8_t* testString, int testStringLength, bool leadingRepeats) {
    if(testStringLength == 0)
        return 0;
    return findNumberOfRepetitions(repeatUnit, repeatUnitLength, 0, repeatUnitLength, testString, testStringLength, 0, testStringLength, leadingRepeats);
}