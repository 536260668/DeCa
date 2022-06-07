//
// Created by lhh on 11/13/21.
//

#include "CommonInfo.h"
#include <cmath>
#include <stdexcept>
#include "VCFConstants.h"

CommonInfo::CommonInfo(std::string & name, double log10PError, std::set<std::string> * filters, std::map<std::string, AttributeValue>* attributes): name(name), attributes(attributes)
{
    setLog10PError(log10PError);
    if(filters != nullptr){
        this->filters = *filters;
    }
}

CommonInfo::~CommonInfo() {
    delete attributes;
}

void CommonInfo::setLog10PError(double log10PError)
{
    if ( log10PError > 0 && log10PError != NO_LOG10_PERROR)
        throw std::invalid_argument("BUG: log10PError cannot be > 0 : ");
    if (std::isinf(this->log10PError))
        throw std::invalid_argument("BUG: log10PError should not be Infinity");
    if ( std::isnan(this->log10PError) )
        throw std::invalid_argument("BUG: log10PError should not be nan");
    this->log10PError = log10PError;
}

bool CommonInfo::hasAttribute(const std::string & key) {
    if(!attributes)
        return false;
    else
        return attributes->find(key) != attributes->end();
}

int CommonInfo::getAttributeAsInt(std::string &key, int defaultValue) {
    return attributes->at(key).getAttributeAsInt();
}

AttributeValue CommonInfo::getAttribute(std::string &key) {
    return attributes->at(key);
}

std::map<std::string, AttributeValue> &CommonInfo::getAttributes(){
    return *attributes;
}

std::set<std::string>* CommonInfo::getFiltersMaybeNull() {
    if(filters.empty())
        return nullptr;
    else
        return &filters;
}

std::set<std::string> &CommonInfo::getFilters()
{
    return filters;
}

bool CommonInfo::filtersWereApplied() {
    return !filters.empty();
}

double CommonInfo::getLog10PError() const {
    return log10PError;
}

std::string & CommonInfo::getName() {
    return name;
}
