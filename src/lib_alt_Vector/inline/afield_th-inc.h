/*!
        @file    afield_th-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#ifndef AFIELD_TH_INC_H
#define AFIELD_TH_INC_H

#ifdef USE_BENCHMARK
#include <omp.h>
#else
#include "lib/ResourceManager/threadManager.h"
#endif

namespace {

//====================================================================
inline void set_thread(int& ith, int& nth)
{
#ifdef USE_BENCHMARK
  nth = omp_get_num_threads();
  ith = omp_get_thread_num();
#else
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();
#endif

}
  
//====================================================================
inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
                           const size_t size)
{
#ifdef USE_BENCHMARK
  nth = omp_get_num_threads();
  ith = omp_get_thread_num();
#else
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();
#endif

  size_t is2 = size * size_t(ith) / nth;
  size_t ns2 = size * size_t(ith + 1) / nth;
  is = int(is2);
  ns = int(ns2);

}

//====================================================================

} // nameless namespace end
#endif
