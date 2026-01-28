/*!
        @file    timer.cpp

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#include "timer.h"

const std::string Timer::class_name = "Timer";

//====================================================================
void Timer::timestamp()
{
  const size_t buf_size = 1024;
  static char  buf[buf_size];

  time_t    current_time;
  struct tm *timep;

  current_time = time(NULL);
  timep        = localtime(&current_time);

  strftime(buf, buf_size, "%Y/%m/%d %H:%M:%S %z", timep);

  vout.general("%s: timestamp: %s\n", class_name.c_str(), buf);
}


//====================================================================
Timer::~Timer()
{
  if (m_report_on_exit) report();
}


//====================================================================
void Timer::start()
{
  struct timeval t_start;

#ifdef USE_RUSAGE
  struct rusage ru;
  int           result = getrusage(RUSAGE_SELF, &ru);
  t_start = ru.ru_utime;
#else
  int result = gettimeofday(&t_start, 0);
#endif

  if (result) {
    vout.general("%s: warning, aquiring system clock failed.\n", class_name.c_str());
    return;
  }

  m_start = (double)t_start.tv_sec + t_start.tv_usec * 1.0e-6;

  is_started = true;
  ++m_counter;
}


//====================================================================
void Timer::stop()
{
  if (!is_started) return;

  struct timeval t_end;

#ifdef USE_RUSAGE
  struct rusage ru;
  int           result = getrusage(RUSAGE_SELF, &ru);
  t_end = ru.ru_utime;
#else
  int result = gettimeofday(&t_end, 0);
#endif

  if (result) {
    vout.general("%s: warning, aquiring system clock failed.\n", class_name.c_str());
    return;
  }

  double m_end = (double)t_end.tv_sec + t_end.tv_usec * 1.0e-6;

  m_elapsed += (m_end - m_start);

  is_started = false;
}


//====================================================================
void Timer::reset()
{
  is_started = false;
  m_elapsed  = double(0);
  m_start    = double(0);
  m_counter  = 0;
}


//====================================================================
double Timer::elapsed_sec() const
{
  return m_elapsed;
}


//====================================================================
double Timer::elapsed_msec() const
{
  return m_elapsed * 1.0e+3;
}


//====================================================================
unsigned long Timer::get_counter() const
{
  return m_counter;
}


//====================================================================
void Timer::report(const Bridge::VerboseLevel vl)
{
  stop();

  unsigned long count   = get_counter();
  double        elapsed = elapsed_sec();
  double        average = count ? elapsed / count : 0.0;

  vout.general(vl, "Elapsed time: %s: total %12.2f sec, count %4d, average %12.2f sec\n", m_id.c_str(), elapsed, count, average);
}


//==========================================================
//==================================================END=====
