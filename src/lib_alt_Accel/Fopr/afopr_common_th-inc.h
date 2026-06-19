/*!
        @file    afopr_common_th-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: aoyama $
        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$
        @version $LastChangedRevision: 1928 $
*/

#ifndef ACCEL_AFOPR_COMMON_TH_INC_INCLUDED
#define ACCEL_AFOPR_COMMON_TH_INC_INCLUDED


namespace {

//====================================================================
// case (a): tasks are equally assgined to all threads
/*
inline void set_threadtask_afopr(int& ith, int& nth, int& is, int& ns,
                           const int size)
{
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();

  is = size * ith / nth;
  ns = size * (ith + 1) / nth;

}
*/
//====================================================================
// case (b): tasks are almost assigned to slave threads.
/*
inline void set_threadtask_afopr(int& ith, int& nth, int& is, int& ns,
                           const int size)
{
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();

  if(nth > 1){
    int offset = size % (nth-1);
    int ntask = size/(nth-1);
    is = offset + ntask * (ith-1);
    if(ith == 0) is = 0;
    ns = offset + ntask * ith;
  }else{
    is = 0;
    ns = size;
  }

}
*/
//====================================================================
// case (c): tasks are assigned only to slave threads.

inline void set_threadtask_afopr(int& ith, int& nth, int& is, int& ns,
                           const int size)
{
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();

  if(nth == 1){
    is = 0;
    ns = size;
  }else{
    is = size * (ith - 1) / (nth-1);
    if(ith == 0) is = 0;
    ns = size * ith / (nth-1);
  }

}


} // nameless namespace end
#endif
//============================================================END=====
