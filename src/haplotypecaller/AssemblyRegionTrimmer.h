//
// Created by 梦想家xixi on 2021/12/14.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_H
#define MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_H

#include "ReadThreadingAssemblerArgumentCollection.h"

class AssemblyRegionTrimmer {
private:
    int usableExtension;
    bool emitReferenceConfidence;
    ReadThreadingAssemblerArgumentCollection* assemblyArgs;
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_H
