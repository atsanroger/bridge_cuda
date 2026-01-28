/*!
        @file    afield_th-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef AFIELD_TH_INC_H
#define AFIELD_TH_INC_H

namespace {

//====================================================================
inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
                             const int size)
{
  nth = ThreadManager_OpenMP::get_num_threads();
  ith = ThreadManager_OpenMP::get_thread_id();

  //  is = size * ith / nth;
  //  ns = size * (ith + 1) / nth;
  size_t is2 = size_t(size) * size_t(ith) / nth;
  size_t ns2 = size_t(size) * size_t(ith + 1) / nth;
  is = int(is2);
  ns = int(ns2);

}

//====================================================================
template <typename REALTYPE>
inline void set_th(REALTYPE *yL, REALTYPE a, int size)
{
  for(int j = 0; j < size; ++j){
    yL[j] = a;
  }

}

//====================================================================
template <typename REALTYPE>
inline void scal_th(REALTYPE *yL, REALTYPE a, int size)
{
  for(int j = 0; j < size; ++j){
    yL[j] *= a;
  }

}

//====================================================================
template <typename REALTYPE>
inline void copy_th(REALTYPE *yL, REALTYPE *xL, int size)
{
  for(int j = 0; j < size; ++j){
    yL[j] = xL[j];
  }

}

//====================================================================
template <typename REALTYPE>
inline void mul_th(REALTYPE *yL, REALTYPE a, REALTYPE *xL, int size)
{
  for(int j = 0; j < size; ++j){
    yL[j] = a * xL[j];
  }

}

//====================================================================
template <typename REALTYPE>
inline void axpy_th(REALTYPE *yL, REALTYPE a, REALTYPE *xL, int size)
{
  for(int j = 0; j < size; ++j){
    yL[j] += a * xL[j];
  }

}

//====================================================================

} // nameless namespace end
#endif
