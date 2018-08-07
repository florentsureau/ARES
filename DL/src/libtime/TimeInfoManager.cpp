
#include "CTimer.h"
#include <limits>
#include <float.h>
#include <thread>
#include <iostream>
#include <sstream>

#include "TimeInfoManager.h"

std::vector<CTimerImpl> 	TimeInfoManager::timersImpl				= std::vector<CTimerImpl>();
std::vector<std::string> 	TimeInfoManager::timersName             = std::vector<std::string>();
std::vector<double> 		TimeInfoManager::timersValue            = std::vector<double>();
std::vector<double> 		TimeInfoManager::timersSum              = std::vector<double>();
std::vector<unsigned int> 	TimeInfoManager::timersCallCounter      = std::vector<unsigned int>();
std::vector<double> 		TimeInfoManager::timersMin              = std::vector<double>();
std::vector<double> 		TimeInfoManager::timersMax              = std::vector<double>();

void TimeInfoManager::_startTimer(unsigned int id)
{
	timersImpl[id].start();
}

void TimeInfoManager::_stopTimer(unsigned int id)
{
	timersImpl[id].stop();
	timersValue[id] = timersImpl[id].seconds();
	timersSum[id] += timersImpl[id].seconds();
	timersCallCounter[id]++;
	if(timersValue[id] < timersMin[id]) timersMin[id] = timersValue[id];
	if(timersValue[id] > timersMax[id]) timersMax[id] = timersValue[id];
}

void TimeInfoManager::_resetTimer(unsigned int id)
{
	timersImpl[id].stop();
	timersValue[id] = 0;
	timersSum[id] = 0;
	timersCallCounter[id] = 0;
	timersMin[id] = DBL_MAX;
	
#ifdef _MSC_VER
	timersMin.push_back( DBL_MAX );
#else
	timersMin.push_back(std::numeric_limits<double>::max() );
#endif

	timersMax[id] = 0;
}

double TimeInfoManager::_getTimerValue(unsigned int id){
	return timersValue[id];
}

double TimeInfoManager::_getTimerSum(unsigned int id){
	return timersSum[id];
}

double TimeInfoManager::_getTimerAverage(unsigned int id){
	return timersSum[id]/timersCallCounter[id];
}

double TimeInfoManager::_getTimerMin(unsigned int id){
	return timersMin[id];
}

double TimeInfoManager::_getTimerMax(unsigned int id){
	return timersMax[id];
}

unsigned int TimeInfoManager::findTimerByName(std::string name)
{
    std::string nameAndThread;
    std::ostringstream oss;
    oss << name << "_" << std::this_thread::get_id();
    nameAndThread = oss.str();
	for(unsigned int i = 0 ; i < timersImpl.size() ; i++)
        if(timersName[i].compare(nameAndThread) == 0)
			return i;

    addTimer(nameAndThread);

	return static_cast<unsigned int>(timersImpl.size()) - 1;
}

void TimeInfoManager::addTimer(std::string name){
	timersImpl.push_back(CTimerImpl());
	timersName.push_back(name);
	timersValue.push_back(0);
	timersSum.push_back(0);
	timersCallCounter.push_back(0);
	
#ifdef _MSC_VER
	timersMin.push_back( DBL_MAX );
#else
	timersMin.push_back(std::numeric_limits<double>::max() );
#endif

	timersMax.push_back(0);
}

void TimeInfoManager::stopTimer(std::string name){
    _stopTimer(findTimerByName(name));
}

void TimeInfoManager::startTimer(std::string name){
    _startTimer(findTimerByName(name));
}

void TimeInfoManager::resetTimer(std::string name){
	_resetTimer(findTimerByName(name));
}

double TimeInfoManager::getTimerValue(std::string name){
	return _getTimerValue(findTimerByName(name));
}

double TimeInfoManager::getTimerSum(std::string name){
	return _getTimerSum(findTimerByName(name));
}

double TimeInfoManager::getTimerAverage(std::string name){
	return _getTimerAverage(findTimerByName(name));
}

double TimeInfoManager::getTimerMin(std::string name){
	return _getTimerMin(findTimerByName(name));
}

double TimeInfoManager::getTimerMax(std::string name){
	return _getTimerMax(findTimerByName(name));
}

void TimeInfoManager::resetAllTimer(){
    for(unsigned int i = 0 ; i < timersImpl.size() ; ++i)
        _resetTimer(i);
}

void TimeInfoManager::display_all()
{

	bool doPercent = false;
	unsigned int idGlobal = 0;
	for(unsigned int i = 0 ; i < timersImpl.size() ; i++)
		if(timersName[i].compare("Global") == 0)
		{
			doPercent = true;
			idGlobal = i;
		}

	for(unsigned int i = 0 ; i < timersImpl.size() ; i++)
	{
        if(0 >= timersCallCounter[i]) {
            continue;
        }
		std::cout << timersName[i] << " : " << std::endl;
		std::cout << "\t" << "Sum = " << _getTimerSum(i) << std::endl;
		std::cout << "\t" << "Average = " << _getTimerAverage(i) << std::endl;
		std::cout << "\t" << "Nb Calls = " << timersCallCounter[i] << std::endl;
		std::cout << "\t" << "Min = " << _getTimerMin(i) << std::endl;
		std::cout << "\t" << "Max = " << _getTimerMax(i) << std::endl;
		if(doPercent)
			std::cout << "\t" << "% = " << std::setprecision(5) << 100 * (_getTimerSum(i)/_getTimerSum(idGlobal) ) << std::endl;
		std::cout << std::endl;
	}

}

void TimeInfoManager::write_all_csv(std::string fileName){

	std::ofstream file;
	file.open(fileName.c_str());

	bool doPercent = false;
	unsigned int idGlobal = 0;
	for(unsigned int i = 0 ; i < timersImpl.size() ; i++)
		if(timersName[i].compare("Global") == 0)
		{
			doPercent = true;
			idGlobal = i;
		}

	file << " ; Sum ; Average ; Nb Calls ; Min ; Max ; % " << std::endl;

	for(unsigned int i = 0 ; i < timersImpl.size() ; i++)
	{
        if(0 >= timersCallCounter[i]) {
            continue;
        }
		file << timersName[i] << " ; " << _getTimerSum(i) << " ; " << _getTimerAverage(i) << " ; ";
		file << timersCallCounter[i] << " ; " << _getTimerMin(i) << " ; " << _getTimerMax(i);
		if(doPercent)
			file << " ; " << _getTimerSum(i)/_getTimerSum(idGlobal) ;
		file << std::endl;
	}

	file.close();

}

void TimeInfoManager::write_all(std::string fileName){

	std::ofstream file;
	file.open(fileName.c_str());
	unsigned int maxLength = 0;

	bool doPercent = false;
	unsigned int idGlobal = 0;
	for(unsigned int i = 0 ; i < timersImpl.size() ; i++)
		if(timersName[i].compare("Global") == 0)
		{
			doPercent = true;
			idGlobal = i;
		}

	for(unsigned int i = 0 ; i < timersName.size() ; i++)
		if(timersName[i].length() > maxLength) maxLength = static_cast<unsigned int>(timersName[i].length());

	for(unsigned int j = 0 ; j < maxLength ; j++)
		file << " ";
	file << "\t| Sum             ";
	file << "\t| Average         ";
	file << "\t| Nb Calls        ";
	file << "\t| Min             ";
	file << "\t| Max             ";
	if(doPercent)
	file << "\t| %               ";
	file << std::endl;

	for(unsigned int i = 0 ; i < timersImpl.size() ; i++)
	{
        if(0 >= timersCallCounter[i]) {
            continue;
        }
		file << timersName[i];
		for(unsigned int j = 0 ; j < (maxLength - timersName[i].length()) ; j++)
			file << " ";
		file << "\t| " << std::setw(15) << std::left << _getTimerSum(i);
		file << "\t| " << std::setw(15) << std::left << _getTimerAverage(i);
		file << "\t| " << std::setw(15) << std::left << timersCallCounter[i];
		file << "\t| " << std::setw(15) << std::left << _getTimerMin(i);
		file << "\t| " << std::setw(15) << std::left << _getTimerMax(i);
		if(doPercent)
		file << "\t| " << std::setw(15) << std::left <<  _getTimerSum(i)/ _getTimerSum(idGlobal);
		file << std::endl;
	}

	file.close();
}

