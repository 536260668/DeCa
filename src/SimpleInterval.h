//
// Created by 梦想家xixi on 2021/10/11.
//

#ifndef MUTECT2CPP_MASTER_SIMPLEINTERVAL_H
#define MUTECT2CPP_MASTER_SIMPLEINTERVAL_H
#include "Mutect2Utils.h"
#include "IntervalUtils.h"
#include <string>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <algorithm>

#include "Locatable.h"

static const char CONTIG_SEPARATOR = ':';
static const char START_END_SEPARATOR = '-';
static const std::string END_OF_CONTIG = "+";

class SimpleInterval : public Locatable
{
private:
    static const long serialVersionUID = 1L;
    int start;
    int end;
    std::string contig;
    bool contiguous(Locatable* other);

public:
    SimpleInterval(const std::string& contig, int start, int end);

    SimpleInterval(SimpleInterval const &simpleInterval);

    /**
     * Makes an interval by parsing the string.
     *
     * @warning this method does not fill in the true contig end values
     * for intervals that reach to the end of their contig,
     * uses {@link Integer#MAX_VALUE} instead.
     *
     * Semantics of start and end are defined in {@link Locatable}.
     * The format is one of:
     *
     * contig           (Represents the whole contig, from position 1 to the {@link Integer#MAX_VALUE})
     *
     * contig:start     (Represents the 1-element range start-start on the given contig)
     *
     * contig:start-end (Represents the range start-end on the given contig)
     *
     * contig:start+    (Represents the prefix of the contig starting at the given start position and ending at {@link Integer#MAX_VALUE})
     *
     * examples (note that _all_ commas in numbers are simply ignored, for human convenience):
     *
     * 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000' or 'chr2:1000000+'
      *
      * @param str non-empty string to be parsed
     */
    SimpleInterval(std::string& str);

    SimpleInterval() {};

    virtual ~SimpleInterval() = default;

    /**
     * Test that these are valid values for constructing a SimpleInterval:
     *    contig cannot be null
     *    start must be >= 1
     *    end must be >= start
     */
    static void validatePositions(const std::string& contig, int start, int end);

    /**
      * Test that these are valid values for constructing a SimpleInterval:
      *    contig cannot be null
      *    start must be >= 1
      *    end must be >= start
      */
    static bool isValid(const std::string& contig, int start, int end);

    /**
     * Parses a number like 100000 or 1,000,000 into an int.
     */
    static int parsePosition(std::string pos);

    bool operator==(const SimpleInterval& interval) const;

    int hashCode() const;

    bool equal(const SimpleInterval& interval) const {return *this == interval;}

    std::string getContig() const override {return contig;}

    int getStart() const override {return start;}

    /**
    * @return the 0-based start position (from the GA4GH spec).
    */
    long getGA4GHStart() const {return start - 1;}

    int getEnd() const override {return end;}

    long getGA4GHEnd() const {return end;}

    int size() const {return end - start + 1;}

    /**
      * Determines whether this interval comes within "margin" of overlapping the provided locatable.
      * This is the same as plain overlaps if margin=0.
      *
      * @param other interval to check
      * @param margin how many bases may be between the two interval for us to still consider them overlapping; must be non-negative
      * @return true if this interval overlaps other, otherwise false
      * @throws IllegalArgumentException if margin is negative
      */
    bool overlapsWithMargin(Locatable *other, int margin) const;

    /**
     * Determines whether this interval overlaps the provided locatable.
     *
     * @param other interval to check
     * @return true if this interval overlaps other, otherwise false
     */
    bool overlaps(Locatable* other) override;

    /**
      * Returns the intersection of the two intervals. The intervals must overlap or IllegalArgumentException will be thrown.
     */
    SimpleInterval intersect(Locatable* other);

    /**
      * Returns a new SimpleInterval that represents the entire span of this and other.  Requires that
      * this and that SimpleInterval are contiguous.
      */
    SimpleInterval mergeWithContiguous(Locatable* other);

    /**
      * Returns a new SimpleInterval that represents the region between the endpoints of this and other.
      *
      * Unlike {@link #mergeWithContiguous}, the two intervals do not need to be contiguous
      *
      * @param other the other interval with which to calculate the span
      * @return a new SimpleInterval that represents the region between the endpoints of this and other.
      */
    SimpleInterval spanWith(Locatable* other);

    /**
      * Returns a new SimpleInterval that represents this interval as expanded by the specified amount in both
      * directions, bounded by the contig start/stop if necessary.
      *
      * @param padding amount to expand this interval
      * @param contigLength length of this interval's contig
      * @return a new SimpleInterval that represents this interval as expanded by the specified amount in both
      *         directions, bounded by the contig start/stop if necessary.
      */
    SimpleInterval* expandWithinContig(int padding, int contigLength);

    //TODO:SimpleInterval* expandWithinContig(int padding, SAMSequenceDictionary sequenceDictionary);

    friend std::ostream & operator<<(std::ostream &os, const SimpleInterval& simpleInterval);
};
#endif //MUTECT2CPP_MASTER_SIMPLEINTERVAL_H