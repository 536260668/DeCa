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


class VariantContext {
private:
    std::string contig;
    hts_pos_t start;
    hts_pos_t stop;
    std::string ID;
    std::shared_ptr<Allele> REF;
    std::shared_ptr<Allele> ALT;
    CommonInfo commonInfo;
    bool fullyDecoded;
    /** A set of the alleles segregating in this context */

    void validateStop();
    bool validate(const std::set<Validation>& validationToPerform);
    static std::vector<std::shared_ptr<Allele>> makeAlleles(std::vector<std::shared_ptr<Allele>> & alleles);
    void validateAlleles();
    void validateGenotypes();
    void determineType();
    void determinePolymorphicType();
    static VariantContextType typeOfBiallelicVariant(const std::shared_ptr<Allele> & ref, const std::shared_ptr<Allele> & allele);

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
    std::vector<std::shared_ptr<Allele>>  & getAlleles();
    static bool hasSymbolicAlleles(const std::vector<std::shared_ptr<Allele>> & alleles);
    std::shared_ptr<Allele> getReference();

    bool hasAllele(const std::shared_ptr<Allele>& allele);
    bool hasAllele(const std::shared_ptr<Allele>& allele, bool ignoreRefState);
    bool hasAllele(const std::shared_ptr<Allele>& allele, bool ignoreRefState, bool considerRefAllele);
    std::vector<std::shared_ptr<Allele>> getAlternateAlleles();
    std::shared_ptr<Allele> getAlternateAllele(int i);
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
                   const std::shared_ptr<std::vector<std::shared_ptr<Allele>>> & alleles,
                   GenoTypesContext* genotypes,
                   double log10PError,
                   std::set<std::string>* filters, std::map<std::string, void*>* attributes,
                   bool fullyDecoded,
                   const std::set<Validation>&  validationToPerform
    );
protected:
    VariantContextType type;
    std::vector<std::shared_ptr<Allele>> alleles;
    GenoTypesContext* genotypes;



};


#endif //MUTECT2CPP_MASTER_VARIANTCONTEXT_H
