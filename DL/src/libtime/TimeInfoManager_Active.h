
#ifndef TIMEINFOMANAGERACTIVE_H_
#define TIMEINFOMANAGERACTIVE_H_

#define Go_Time(name)       TimeInfoManager::startTimer(name);
#define Og_Time(name)       TimeInfoManager::stopTimer(name);
#define Time_Display_All    TimeInfoManager::display_all();
#define Time_Reset_All      TimeInfoManager::resetAllTimer();

#endif /* TIMEINFOMANAGERACTIVE_H_ */
