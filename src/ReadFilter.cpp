//
// Created by lhh on 10/19/21.
//

#include "ReadFilter.h"
#include "Mutect2Engine.h"

bool ReadFilter::ReadLengthTest(bam1_t *read)
{
    if(read->core.l_qseq >= Mutect2Engine::MIN_READ_LENGTH)
        return true;
    else
        return false;
}