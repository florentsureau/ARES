
#ifndef TIMEINFOMANAGER_H_
#define TIMEINFOMANAGER_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
class CTimerImpl;


class TimeInfoManager {
private:
    static std::vector<CTimerImpl> timersImpl;
    static std::vector<std::string> timersName;
    static std::vector<double> timersValue;
    static std::vector<double> timersSum;
    static std::vector<unsigned int> timersCallCounter;
    static std::vector<double> timersMin;
    static std::vector<double> timersMax;


	static void _startTimer(unsigned int id);

	static void _stopTimer(unsigned int id);

	static void _resetTimer(unsigned int id);

	static double _getTimerValue(unsigned int id);

	static double _getTimerSum(unsigned int id);

	static double _getTimerAverage(unsigned int id);

	static double _getTimerMin(unsigned int id);

	static double _getTimerMax(unsigned int id);

    static unsigned int findTimerByName(std::string name);

public:

    static void resetAllTimer();

	static void display_all();

	static void write_all_csv(std::string fileName);

	static void write_all(std::string fileName);

    static void addTimer(std::string name);

    static void stopTimer(std::string name);

    static void startTimer(std::string name);

    static void resetTimer(std::string name);

    static double getTimerValue(std::string name);

    static double getTimerSum(std::string name);

    static double getTimerAverage(std::string name);

    static double getTimerMin(std::string name);

    static double getTimerMax(std::string name);

};


#ifdef ENABLE_TIMEINFOMANAGER
    #include "TimeInfoManager_Active.h"
#else
    #include "TimeInfoManager_Passive.h"
#endif /* ENABLE_TIMEINFOMANAGER */

#endif /* TIMEINFOMANAGER_H_ */
