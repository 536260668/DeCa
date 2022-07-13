//
// Created by hlf on 7/8/22.
//

#include "PairHMMConcurrentControl.h"

std::priority_queue<std::shared_ptr<LikelihoodsTask>> PairHMMConcurrentControl::pairHMMTaskQueue;
bool PairHMMConcurrentControl::startPairHMMConcurrentMode;
std::mutex PairHMMConcurrentControl::pairHMMMutex;

void PairHMMConcurrentControl::initial() {
	pairHMMTaskQueue = std::priority_queue<std::shared_ptr<LikelihoodsTask>>();
	startPairHMMConcurrentMode = false;
}

LikelihoodsTask::LikelihoodsTask(std::vector<testcase> *taskTestcases, std::vector<double> *taskLikelihoodArray,
                                 unsigned long index, unsigned long testcasesSize)
		: taskTestcases(taskTestcases), taskLikelihoodArray(taskLikelihoodArray), index(index),
		  count(index), testcasesSize(testcasesSize) {}
