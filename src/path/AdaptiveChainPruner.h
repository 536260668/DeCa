//
// Created by 梦想家xixi on 2021/10/29.
//

#ifndef MUTECT2CPP_MASTER_ADAPTIVECHAINPRUNER_H
#define MUTECT2CPP_MASTER_ADAPTIVECHAINPRUNER_H

#include "Mutect2Utils.h"
#include <vector>
#include <set>
#include "./Path.h"
#include "ChainPruner.h"

template<class V, class E>
class AdaptiveChainPruner : public ChainPruner<V, E>{
private:
    double initialErrorProbability;
    double logOddsThreshold;
    int maxUnprunedVariants;

public:
    AdaptiveChainPruner(const double initialErrorProbability, const double logOddsThreshold, const int maxUnprunedVariants) : ChainPruner<V, E>(),
    initialErrorProbability(initialErrorProbability), logOddsThreshold(logOddsThreshold), maxUnprunedVariants(maxUnprunedVariants) {
        Mutect2Utils::validateArg(initialErrorProbability > 0, "Must have positive error probability");
    }

protected:
    std::unordered_set<Path<V,E>*> chainsToRemove(std::vector<Path<V,E>*> chains) {
        if(chains.empty()){
            std::unordered_set<Path<V,E>*> result;
            return result;
        }

        std::shared_ptr<DirectedSpecifics<V, E>> graph = chains[0]->getGraph();
        std::unordered_set<Path<V,E>*> probableErrorChains = likelyErrorChains(chains, graph, 0.001);
        int errorCount = 0;
        int totalBases = 0;
        typename std::vector<Path<V,E>*>::iterator viter;
        typename std::unordered_set<Path<V,E>*>::iterator siter;
        for(siter = probableErrorChains.begin(); siter != probableErrorChains.end(); siter++) {
            errorCount += (*siter)->getLastEdge()->getMultiplicity();
        }
        for(viter = chains.begin(); viter != chains.end(); viter++) {
            for(typename std::vector<std::shared_ptr<E>>::iterator iter = (*viter)->getEdges().begin(); iter != (*viter)->getEdges().end(); iter++) {
                totalBases += (*iter)->getMultiplicity();
            }
        }
        double errorRate = (double) errorCount / totalBases;
        return likelyErrorChains(chains, graph, errorRate);
    }

private:
    std::unordered_set<Path<V,E>*> likelyErrorChains(std::vector<Path<V,E>*> & chains, std::shared_ptr<DirectedSpecifics<V, E>> graph, double errorRate) {
        std::map<Path<V,E>*, double> chainLogOddsmap;
        typename std::vector<Path<V,E>*>::iterator viter;
        for(viter = chains.begin(); viter != chains.end(); viter++) {
            chainLogOddsmap.insert(std::pair<Path<V,E>* , double>(*viter, chainLogOdds(*viter, graph, errorRate)));
        }
        std::unordered_set<Path<V,E>*> result;
        typename std::map<Path<V,E>*, double>::iterator miter;
        for(miter = chainLogOddsmap.begin(); miter != chainLogOddsmap.end(); miter++) {
            if(miter->second < 2.302585092994046) {
                result.insert(miter->first);
            }
        }
        std::vector<Path<V,E>*> newchains;
        for(viter = chains.begin(); viter != chains.end(); viter++) {
            if(isChainPossibleVariant(*viter, graph))
                newchains.template emplace_back(*viter);
        }
        std::sort(newchains.begin(), newchains.end(), [chainLogOddsmap](Path<V,E>* a, Path<V,E>* b)->bool{return chainLogOddsmap.at(a) > chainLogOddsmap.at(b);});
        std::sort(newchains.begin(), newchains.end(), [](Path<V,E>* a, Path<V,E>* b)->bool{return a->length() < b->length();});
        if(newchains.size() <= 100)
            return result;
        else {
            for(viter = newchains.begin() + 100; viter != newchains.end(); viter++) {
                result.insert(*viter);
            }
            return result;
        }
    }

    double chainLogOdds(Path<V,E>* chain, std::shared_ptr<DirectedSpecifics<V, E>> graph, double errorRate) {
        typename std::unordered_set<std::shared_ptr<E>>::iterator eiter;
        typename std::vector<std::shared_ptr<E>>::iterator viter;
        for(viter = chain->getEdges().begin(); viter != chain->getEdges().end(); viter++) {
            if((*viter)->getIsRef())
                return POSITIVE_INFINITY;
        }
        int leftTotalMultiplicity = 0;
        int rightTotalMultiplicity = 0;
        std::unordered_set<std::shared_ptr<E>> outgoing = graph->outgoingEdgesOf(chain->getFirstVertex());
        std::unordered_set<std::shared_ptr<E>> incoming = graph->incomingEdgesOf(chain->getLastVertex());
        for(eiter = outgoing.begin(); eiter != outgoing.end(); eiter++) {
            leftTotalMultiplicity += (*eiter)->getMultiplicity();
        }
        for(eiter= incoming.begin(); eiter != incoming.end(); eiter++) {
            rightTotalMultiplicity += (*eiter)->getMultiplicity();
        }
        int leftMultiplicity = (chain->getEdges()[0])->getMultiplicity();
        int rightMultiplicity = chain->getLastEdge()->getMultiplicity();

        double leftLogOdds = graph->isSource(chain->getFirstVertex()) ? 0.0 : Mutect2Utils::logLikelihoodRatio(leftTotalMultiplicity - leftMultiplicity, leftMultiplicity, errorRate);
        double rightLogOdds = graph->isSink(chain->getLastVertex()) ? 0.0 : Mutect2Utils::logLikelihoodRatio(rightTotalMultiplicity - rightMultiplicity, rightMultiplicity, errorRate);

        return std::max(leftLogOdds, rightLogOdds);
    }

    bool isChainPossibleVariant(Path<V,E>* chain, std::shared_ptr<DirectedSpecifics<V, E>> graph) {
        typename std::unordered_set<std::shared_ptr<E>>::iterator eiter;
        typename std::vector<std::shared_ptr<E>>::iterator viter;
        for(viter = chain->getEdges().begin(); viter != chain->getEdges().end(); viter++) {
            if((*viter)->getIsRef())
                return POSITIVE_INFINITY;
        }
        int leftTotalMultiplicity = 0;
        int rightTotalMultiplicity = 0;
        std::unordered_set<std::shared_ptr<E>> outgoing = graph->outgoingEdgesOf(chain->getFirstVertex());
        std::unordered_set<std::shared_ptr<E>> incoming = graph->outgoingEdgesOf(chain->getFirstVertex());
        for(eiter = outgoing.begin(); eiter != outgoing.end(); eiter++) {
            leftTotalMultiplicity += (*eiter)->getMultiplicity();
        }
        for(eiter= outgoing.begin(); eiter != outgoing.end(); eiter++) {
            rightTotalMultiplicity += (*eiter)->getMultiplicity();
        }
        int leftMultiplicity = (chain->getEdges()[0])->getMultiplicity();
        int rightMultiplicity = chain->getLastEdge()->getMultiplicity();

        return leftMultiplicity <= leftTotalMultiplicity / 2 || rightMultiplicity <= rightTotalMultiplicity / 2;
    }
};


#endif //MUTECT2CPP_MASTER_ADAPTIVECHAINPRUNER_H
