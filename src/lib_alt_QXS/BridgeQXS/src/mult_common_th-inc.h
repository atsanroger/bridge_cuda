/*!
      @file    mult_common_th-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef QXS_MULT_COMMON_TH_INCLUDED
#define QXS_MULT_COMMON_TH_INCLUDED

#ifdef USE_BENCHMARK
#include <omp.h>
#else
#include "lib/ResourceManager/threadManager.h"
#endif

namespace {
//====================================================================
// case (a): tasks are equally assgined to all threads

  inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
                             const int size)
  {
    nth = ThreadManager::get_num_threads();
    ith = ThreadManager::get_thread_id();

    is = size * ith / nth;
    ns = size * (ith + 1) / nth;
  }


//====================================================================
// case (b): tasks are almost assigned to slave threads.

/*
inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
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

/*
inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
                           const int size)
{
#ifdef USE_BENCHMARK
  nth = omp_get_num_threads();
  ith = omp_get_thread_num();
#else
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();
#endif

  if(nth == 1){
    is = 0;
    ns = size;
  }else{
    is = size * (ith - 1) / (nth-1);
    if(ith == 0) is = 0;
    ns = size * ith / (nth-1);
  }

}
*/
  inline void set_threadtask_without_master(int& ith, int& nth, int& is, int& ns,
                                            const int size)
  {
#ifdef USE_BENCHMARK
    nth = omp_get_num_threads();
    ith = omp_get_thread_num();
#else
    nth = ThreadManager::get_num_threads();
    ith = ThreadManager::get_thread_id();
#endif

    if (nth == 1) {
      is = 0;
      ns = size;
    } else {
      if (ith == 0) {
        is = 0;
        ns = 0;
      } else {
        int ith_tmp = ith - 1;
        int nth_tmp = nth - 1;
        if (nth_tmp == 1) {
          is = 0;
          ns = size;
        } else {
          is = size * (ith_tmp - 1) / (nth_tmp - 1);
          if (ith_tmp == 1) is = 0;
          ns = size * ith_tmp / (nth_tmp - 1);
        }
      }
    }
  }


//====================================================================
} // nameless namespace end

#endif
//============================================================END=====
