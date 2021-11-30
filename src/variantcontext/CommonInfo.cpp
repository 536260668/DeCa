//
// Created by lhh on 11/13/21.
//

#include "CommonInfo.h"
#include <cmath>

CommonInfo::CommonInfo(std::string name, double log10PError, std::set<std::string> filters): name(name), filters(filters)
{
    setLog10PError(log10PError);

}

void CommonInfo::setLog10PError(double log10PError)
{
    if ( log10PError > 0 && log10PError != NO_LOG10_PERROR)
        throw "BUG: log10PError cannot be > 0 : ";
    if (std::isinf(this->log10PError))
        throw "BUG: log10PError should not be Infinity";
    if ( std::isnan(this->log10PError) )
        throw "BUG: log10PError should not be NaN";
    this->log10PError = log10PError;
}
