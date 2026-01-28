/*!
        @file    afield_th-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef ACCEL_AFIELD_TH_INC_H
#define ACCEL_AFIELD_TH_INC_H

#include <cuda_runtime.h>

namespace {

//====================================================================
inline void set_thread(int& ith, int& nth)
{
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();
}

//====================================================================
inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
                           const int size)
{
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();

  size_t is2 = size_t(size) * size_t(ith) / nth;
  size_t ns2 = size_t(size) * size_t(ith + 1) / nth;
  is = int(is2);
  ns = int(ns2);

}

//====================================================================

} // nameless namespace end
#endif
