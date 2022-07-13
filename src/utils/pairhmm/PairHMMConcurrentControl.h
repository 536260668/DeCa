//
// Created by hlf on 7/8/22.
//

#ifndef MUTECT2CPP_MASTER_PAIRHMMCONCURRENTCONTROL_H
#define MUTECT2CPP_MASTER_PAIRHMMCONCURRENTCONTROL_H

#include "pairhmm_common.h"
#include <vector>
#include <atomic>
#include <mutex>
#include <queue>
#include <condition_variable>

struct LikelihoodsTask {
	std::vector<testcase> *taskTestcases = nullptr;
	std::vector<double> *taskLikelihoodArray = nullptr;
	std::atomic<unsigned long> count = 0;
	std::atomic<unsigned long> index = 0;
	unsigned long testcasesSize = 0;

	LikelihoodsTask(std::vector<testcase> *taskTestcases, std::vector<double> *taskLikelihoodArray, unsigned long index,
	                unsigned long testcasesSize);

	friend bool operator < (const LikelihoodsTask& a, const LikelihoodsTask & b){
		return a.testcasesSize < b.testcasesSize;
	}
};

class PairHMMConcurrentControl {
public:
	static std::mutex pairHMMMutex;
	static std::priority_queue<std::shared_ptr<LikelihoodsTask>> pairHMMTaskQueue;
	static bool startPairHMMConcurrentMode;

	static void initial();
};


#endif //MUTECT2CPP_MASTER_PAIRHMMCONCURRENTCONTROL_H
