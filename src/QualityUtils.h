/**
 * QualityUtils is a static class with some utility methods for manipulating quality scores
 */

#ifndef QUALITY_UTILS_H
#define QUALITY_UTILS_H


class QualityUtils{
public:
    const static char MIN_USABLE_Q_SCORE = 6;
    const static int MAPPING_QUALITY_UNAVALIABLE = 255;

    /**
     * bams containing quals above this value are extremely suspicious and we should warn the user
     */
    const static char MAX_REASONABLE_Q_SCORE = 60;

    const static char MAX_SAM_QUAL_SCORE = 93;

    static double qualToErrorProb(int qual);

    /**
     * Convert a probability of being wrong to a phred-scaled quality score (0.01 => 20).
     *
     * Note, this function caps the resulting quality score by the public static value MAX_SAM_QUAL_SCORE
     * and by 1 at the low-end.
     *
     * @param errorRate a probability (0.0-1.0) of being wrong (i.e., 0.01 is 1% change of being wrong)
     * @return a quality score (0-MAX_SAM_QUAL_SCORE)
     */
    static char errorProbToQual(double errorRate);

    /**
     * Convert a probability of being wrong to a phred-scaled quality score (0.01 => 20).
     *
     * Note, this function caps the resulting quality score by the public static value MIN_REASONABLE_ERROR
     * and by 1 at the low-end.
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param errorRate a probability (0.0-1.0) of being wrong (i.e., 0.01 is 1% change of being wrong)
     * @return a quality score (0-maxQual)
     */
    static char errorProbToQual(double errorRate, char maxQual);

    static char boundQual(int qual, char maxQual);
};

#endif