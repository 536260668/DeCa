//
// Created by lhh on 4/26/22.
//

#ifndef MUTECT2CPP_MASTER_GATKVARIANTCONTEXTUTILS_H
#define MUTECT2CPP_MASTER_GATKVARIANTCONTEXTUTILS_H


#include <cstdint>

class GATKVariantContextUtils {
public:
    /**
     * Finds number of repetitions a string consists of.
     * Same as {@link #findNumberOfRepetitions} but operates on subarrays of a bigger array to save on copying.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnitFull             Non-empty substring represented by byte array
     * @param repeatUnitFullLength       the total length of repeatUnitFull
     * @param offsetInRepeatUnitFull     the offset in repeatUnitFull from which to read the repeat unit
     * @param repeatUnitLength           length of the repeat unit
     * @param testStringFull             string to test (represented by byte array), may be empty
     * @param testStringFullLength       the total length of testStringFull
     * @param offsetInTestStringFull     the offset in offsetInRepeatUnitFull from which to read the test string
     * @param testStringLength           length of the test string
     * @param leadingRepeats         Look for leading (at the beginning of string) or trailing (at end of string) repetitions
     * For example:
     *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
     *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
     *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
     * @return  Number of repetitions (0 if testString is not a concatenation of n repeatUnit's, including the case of empty testString)
     */
    static int findNumberOfRepetitions(uint8_t* repeatUnitFull, int repeatUnitFullLength, int offsetInRepeatUnitFull, int repeatUnitLength, uint8_t* testStringFull, int testStringFullLength, int offsetInTestStringFull, int testStringLength, bool leadingRepeats);

    static int findNumberOfRepetitions(uint8_t* repeatUnit, int repeatUnitLength, uint8_t* testString, int testStringLength, bool leadingRepeats);
};


#endif //MUTECT2CPP_MASTER_GATKVARIANTCONTEXTUTILS_H
