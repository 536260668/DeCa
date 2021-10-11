//
// Created by 梦想家xixi on 2021/10/11.
//

#include "SimpleInterval.h"

int main(int argc, char **argv)
{
    SimpleInterval* interval = new SimpleInterval(1, 10, "123");
    interval->printfInterval();
    delete interval;
    return 0;
}
