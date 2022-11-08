//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_BETADISTRIBUTIONSHAPE_H
#define MUTECT2CPP_MASTER_BETADISTRIBUTIONSHAPE_H


class BetaDistributionShape {
private:
    double alpha;
    double beta;

public:
    static const BetaDistributionShape FLAT_BETA;

    BetaDistributionShape(double alpha, double beta);
    double getAlpha() const;
    double getBeta() const;
};


#endif //MUTECT2CPP_MASTER_BETADISTRIBUTIONSHAPE_H
