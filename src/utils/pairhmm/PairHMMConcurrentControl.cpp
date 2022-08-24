//
// Created by hlf on 7/8/22.
//

#include "PairHMMConcurrentControl.h"

std::priority_queue<std::shared_ptr<LikelihoodsTask>> PairHMMConcurrentControl::pairHMMTaskQueue;
bool PairHMMConcurrentControl::startPairHMMConcurrentMode;
std::mutex PairHMMConcurrentControl::pairHMMMutex;
//std::atomic<unsigned long long> PairHMMConcurrentControl::all_cases;
//std::atomic<unsigned long long> PairHMMConcurrentControl::unique_cases;
//std::atomic<unsigned long long> PairHMMConcurrentControl::all_reads;
//std::atomic<unsigned long long> PairHMMConcurrentControl::unique_reads;

void PairHMMConcurrentControl::initial() {
	pairHMMTaskQueue = std::priority_queue<std::shared_ptr<LikelihoodsTask>>();
	startPairHMMConcurrentMode = false;
//	all_cases = 0;
//	unique_cases = 0;
//	all_reads = 0;
//	unique_reads = 0;
}

LikelihoodsTask::LikelihoodsTask(std::vector<testcase> *taskTestcases, std::vector<double> *taskLikelihoodArray,
                                 unsigned long index, unsigned long testcasesSize)
		: taskTestcases(taskTestcases), taskLikelihoodArray(taskLikelihoodArray), index(index),
		  count(index), testcasesSize(testcasesSize) {}
