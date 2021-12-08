//
// Created by lhh on 11/11/21.
//

#ifndef MUTECT2CPP_MASTER_VARIANTCONTEXT_H
#define MUTECT2CPP_MASTER_VARIANTCONTEXT_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include "htslib/sam.h"
#include "haplotype/Allele.h"
#include "VCFConstants.h"
#include "GenoTypesContext.h"
#include "CommonInfo.h"

class GenoTypesContext;

enum VariantContextType{VariantContext_NO_VARIATION,
    VariantContext_SNP,
    VariantContext_MNP,
    VariantContext_INDEL,
    VariantContext_SYMBOLIC,
    VariantContext_MIXED,
    VariantContext_NULL};

enum Validation{ALLELES,
    GENOTYPES};

//TODO: add GenoTypesContext class 2021.11.11
class VariantContext {
private:
    std::string contig;
    hts_pos_t start;
    hts_pos_t stop;
    std::string ID;
    Allele* REF;
    Allele* ALT;
    CommonInfo commonInfo;
    bool fullyDecoded;
    /** A set of the alleles segregating in this context */

    void validateStop();
    bool validate(const std::set<Validation>& validationToPerform);
    static std::vector<Allele*> makeAlleles(std::vector<Allele*> & alleles);
    void validateAlleles();
    void validateGenotypes();
    void determineType();
    void determinePolymorphicType();
    static VariantContextType typeOfBiallelicVariant(Allele* ref, Allele* allele);

public:
    /**
     * the actual constructor.  Private access only
     *
     * @param source          source
     * @param contig          the contig
     * @param start           the start base (one based)
     * @param stop            the stop reference base (one based)
     * @param alleles         alleles
     * @param genotypes       genotypes map
     * @param log10PError  qual
     * @param filters         filters: use null for unfiltered and empty set for passes filters
     * @param attributes      attributes
     * @param validationToPerform     set of validation steps to take
     */

    bool hasAttribute(std::string &key);

    int getAttributeAsInt(std::string &key, int defaultValue);

    int getEnd() const;
    int getStart();
    bool isBiallelic();
    int getNAlleles();
    bool isSNP();
    bool isSimpleDeletion();
    bool isSimpleInsertion();
    bool isSimpleIndel();
    VariantContextType getType();

    bool hasSymbolicAlleles();
    std::vector<Allele*>  & getAlleles();
    static bool hasSymbolicAlleles(std::vector<Allele*> & alleles);
    Allele* getReference();

    bool hasAllele(Allele* allele);
    bool hasAllele(Allele* allele, bool ignoreRefState);
    bool hasAllele(Allele* allele, bool ignoreRefState, bool considerRefAllele);
    std::vector<Allele*> getAlternateAlleles();
    Allele* getAlternateAllele(int i);
    std::map<std::string, void*> & getAttributes();
    std::string & getContig();
    std::set<std::string>  * getFiltersMaybeNull();
    GenoTypesContext* getGenotypes();
    std::string & getID();
    double getLog10PError();
    std::string  & getSource();
    bool isFullyDecoded();

    VariantContext(std::string &source,
                   std::string &ID,
                   std::string &contig,
                   long start,
                   long stop,
                   std::vector<Allele*> *alleles,
                   GenoTypesContext* genotypes,
                   double log10PError,
                   std::set<std::string>* filters, std::map<std::string, void*>* attributes,
                   bool fullyDecoded,
                   std::set<Validation> & validationToPerform
    );
protected:
    VariantContextType type;
    std::vector<Allele*>  alleles;
    GenoTypesContext* genotypes;



};


#endif //MUTECT2CPP_MASTER_VARIANTCONTEXT_H
