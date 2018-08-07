/******************************************************************************
 **                   Copyright (C) 2017 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Florent Sureau
 **
 **    Date:  01-02/2017
 **
 **    File:  DLException.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Helper routines for exception messages
 **    -----------
 ******************************************************************************/
 
#ifndef DLEXCEPTION_H_
#define DLEXCEPTION_H_
#include <vector>
#include <string.h>
#include <exception>
#include <sstream>
#include <typeinfo>
#include <string>

static std::string createExceptionMsg(const char* func, const char *file, int line) {
    std::stringstream strstr;
    std::string exMessage(func);
    exMessage.append("- in file ");
    strstr<<file;
    exMessage.append(strstr.str());
    strstr.clear();
    exMessage.append(", line ");
    strstr<<line;
    exMessage.append(strstr.str());
    exMessage.append(":");
    strstr<<std::endl;
    return exMessage;
}


#define EXCEPT_WHERE() createExceptionMsg(__PRETTY_FUNCTION__,__FILE__,__LINE__)

#endif //DLEXCEPTION_H_
