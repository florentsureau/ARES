#ifndef _TIMER_H_
#define _TIMER_H_

//--------------------------------------------- windows performance counter ----
#if defined(_WIN32) && (defined(__MINGW32__) || defined(_MSC_VER) || defined(__INTEL_COMPILER))

#include <windows.h>
#include <stdio.h>

class CTimerImpl
{
protected:
   LARGE_INTEGER freq;
   LARGE_INTEGER count;
   LARGE_INTEGER startTime;

public:
   CTimerImpl()
   {
       if (QueryPerformanceFrequency(&freq)==FALSE)
       {
          printf("Performance counter unavailable\n");
          throw;
       }
       reset();
   }

   void reset()
   {
       memset(&count,    0, sizeof(count));
       memset(&startTime,0, sizeof(count));
   }
   void start()
   {
       reset();
       QueryPerformanceCounter(&startTime);
   }
   void stop()
   {
       LARGE_INTEGER stopTime;

       QueryPerformanceCounter(&stopTime);
       count.QuadPart += stopTime.QuadPart - startTime.QuadPart;
   }
   void resume()
   {
       QueryPerformanceCounter(&startTime);
   }
   double seconds() const
   {
       return (double)count.QuadPart/(double)freq.QuadPart;
   }
};


//-------------------------------------------------------------- posix time ----
#elif defined(__GNUC__) && defined(__POSIX__)

#include <ctime>

class CTimerImpl
{
protected:
   clockid_t id;
   double    elapsedTime;
   timespec  startTime;

public:
    CTimerImpl(clockid_t iId) :
        id(iId)
    {
        reset();
    }

    void reset()
    {
        elapsedTime = 0.0;
    }
    void start()
    {
        reset();
        clock_gettime(id, &startTime);
    }
    void stop()
    {
        timespec stop;
        clock_gettime(id, &stopTime);
        elapsedTime += stopTime.tv_sec - startTime.tv_sec;
        elapsedTime += double(stopTime.tv_nsec-startTime.tv_nsec) * 1e-9;
    }
    void resume()
    {
        clock_gettime(id, &startTime);
    }
    double seconds() const
    {
        return elapsedTime;
    }
};

//------------------------------------------------------------ gettimeofday ----
#elif defined(__GNUC__) || (defined(__INTEL_COMPILER) && !defined(_WIN32))

#include <sys/time.h>
#include <unistd.h>

class CTimerImpl
{
protected:
  struct timeval  startTime;
  struct timeval  stopTime;
  struct timezone zone;

  double elapsedTime;

public:
    CTimerImpl()
    {
        reset();
    }

    void reset()
    {
        elapsedTime = 0.0;
    }
    void start()
    {
        reset();
        gettimeofday(&startTime, &zone);
    }
    void stop()
    {
        gettimeofday(&stopTime, &zone);

        elapsedTime += double(stopTime.tv_sec  - startTime.tv_sec);
        elapsedTime += double(stopTime.tv_usec - startTime.tv_usec) * 1e-6;
    }
    void resume()
    {
        gettimeofday(&startTime, &zone);
    }
    double seconds() const
    {
        return elapsedTime;
    }
};


//------------------------------------------------- standard implementation ----
#else

#include <ctime>

static const unsigned long clockticks = CLOCKS_PER_SEC;

class CTimerImpl
{
protected:
   unsigned long freq;
   unsigned long count;
   unsigned long startTime;

public:
    TimerImplStd() :
        freq(clockticks)
    {
        reset();
    }

    void reset()
    {
        count = 0;
    }
    void start()
    {
        reset();
        startTime = clock();
    }
    void stop()
    {
        unsigned long stopTime = clock();
        count += stopTime - startTime;
    }
    void resume()
    {
        startTime = clock();
    }
    double seconds() const
    {
        return (double)count / (double)freq;
    }
};

#endif //Timer platform specific implementations





//--------------------------------------------------------------- CTimer API ----

/**
 A timer class based almost entirely on the one from OpenMesh
 */
class CTimer
{
protected:
    enum State
    {
        Invalid = -1,
        Stopped =  0,
        Running =  1
    };

    CTimerImpl impl;
    State state;

public:
    CTimer();

    bool is_valid() const;
    bool is_stopped() const;

    void reset();
    void start();
    void stop();
    void resume();

    //@{ (if the timer is in the state 'Stopped')
    /// Returns measured time in seconds
    double seconds() const;
    /// Returns measured time in hundredth seconds
    double hseconds() const;
    /// Returns measured time in milli-seconds
    double mseconds() const;
    /// Returns measured time in micro seconds
    double useconds() const;
    //@}

    //@{ Compare timer values
    bool operator  <(const CTimer& t2) const;
    bool operator  >(const CTimer& t2) const;
    bool operator ==(const CTimer& t2) const;
    bool operator <=(const CTimer& t2) const;
    bool operator >=(const CTimer& t2) const;
    //@}
};


//------------------------------------------------ Timer API implementation ----

#include <cassert>

inline
CTimer::
CTimer() :
#if defined(_WIN32) && defined(_MSC_VER)
    impl(),
#elif defined(__GNUC__) && defined(__POSIX__)
// CLOCK_REALTIME
// CLOCK_MONOTONIC     - ?
// CLOCK_REALTIME_HR   - RTlinux
// CLOCK_MONOTONIC_HR  - ?
  #if defined(CLOCK_REALTIME_HR)
    impl(CLOCK_REALTIME_HR),
  #else
    impl(CLOCK_REALTIME),
  #endif
#elif defined(__GNUC__) || (defined(__INTEL_COMPILER) && !defined(_WIN32))
    impl(),
#else
    impl(),
#endif
    state(Stopped)
{}

inline
bool CTimer::
is_valid() const
{
    return state!=Invalid;
}

inline
bool CTimer::
is_stopped() const
{
    return state==Stopped;
}

inline
void CTimer::
reset()
{
    impl.reset();
    state = Stopped;
}

inline
void CTimer::
start()
{
    impl.start();
    state = Running;
}

inline
void CTimer::
stop()
{
    if (state == Running)
    {
        impl.stop();
        state = Stopped;
    }
}

inline
void CTimer::
resume()
{
    if (state == Stopped)
    {
        impl.resume();
        state = Running;
    }
}

inline
double CTimer::
seconds() const
{
    return state==Stopped ? impl.seconds() : 0.0;
}

inline
double CTimer::
hseconds() const
{
    return seconds() * 1e2;
}

inline
double CTimer::
mseconds() const
{
    return seconds() * 1e3;
}

inline
double CTimer::
useconds() const
{
    return seconds() * 1e6;
}

inline
bool CTimer::
operator <(const CTimer& t2) const
{
    assert(is_stopped() && t2.is_stopped());
    return seconds()<t2.seconds();
}

inline
bool CTimer::
operator >(const CTimer& t2) const
{
    assert(is_stopped() && t2.is_stopped());
    return seconds()>t2.seconds();
}

inline
bool CTimer::
operator == (const CTimer& t2) const
{
    assert(is_stopped() && t2.is_stopped());
    return (int)mseconds() == (int)t2.mseconds();
}

inline
bool CTimer::
operator <=(const CTimer& t2) const
{
    assert(is_stopped() && t2.is_stopped());
    return seconds()<=t2.seconds();
}

inline
bool CTimer::
operator >=(const CTimer& t2) const
{
    assert(is_stopped() && t2.is_stopped());
    return seconds()>=t2.seconds();
}

#endif //_TIMER_H_
