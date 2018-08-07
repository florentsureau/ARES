/******************************************************************************
 **                   Copyright (C) 2016 by CEA/Safran
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Arnaud Woiselle
 **
 **    Date:  2016/12
 **
 **    File:  myexit.hpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION  Compliance class for exit(int) function with mex
 **    -----------
 **
 **
 ******************************************************************************/

#ifndef MYEXIT_H_
#define MYEXIT_H_

#include <iostream>
class MyExit {
public:
    static void (*exitPtr)(int);
    static void exit(int errID) {
        exitPtr(errID);
    }
};

#endif
