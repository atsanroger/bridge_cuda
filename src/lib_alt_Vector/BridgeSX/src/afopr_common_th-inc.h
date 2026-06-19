/*!
      @file    afopr_common_th-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#ifndef BRIDGEVEC_AFOPR_COMMON_TH_INCLUDED
#define BRIDGEVEC_AFOPR_COMMON_TH_INCLUDED


#include "lib/ResourceManager/threadManager.h"

namespace {

//====================================================================
// case (a): tasks are equally assgined to all threads
/*
inline void set_threadtask_mult(int& ith, int& nth, int& is, int& ns,
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
inline void set_threadtask_mult(int& ith, int& nth, int& is, int& ns,
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

inline void set_threadtask_mult(int& ith, int& nth, int& is, int& ns,
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

//====================================================================
template <typename REALTYPE>
inline void mult_gm5_dirac(REALTYPE *yL, REALTYPE *xL, int Nc)
{
  for(int ic = 0; ic < Nc; ++ic){
    for(int j = 0; j < VLEN; ++j){
      yL[j + VLEN*(ic + Nc* 0)] = xL[j + VLEN*(ic + Nc* 2)];
      yL[j + VLEN*(ic + Nc* 1)] = xL[j + VLEN*(ic + Nc* 3)];
      yL[j + VLEN*(ic + Nc* 2)] = xL[j + VLEN*(ic + Nc* 0)];
      yL[j + VLEN*(ic + Nc* 3)] = xL[j + VLEN*(ic + Nc* 1)];
    }
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_gm5_chiral(REALTYPE *yL, REALTYPE *xL, int Nc)
{
  for(int ic = 0; ic < Nc; ++ic){
    for(int j = 0; j < VLEN; ++j){
      yL[j + VLEN*(ic + Nc* 0)] =   xL[j + VLEN*(ic + Nc* 0)];
      yL[j + VLEN*(ic + Nc* 1)] =   xL[j + VLEN*(ic + Nc* 1)];
      yL[j + VLEN*(ic + Nc* 2)] = - xL[j + VLEN*(ic + Nc* 2)];
      yL[j + VLEN*(ic + Nc* 3)] = - xL[j + VLEN*(ic + Nc* 3)];
    }
  }

}

//====================================================================
template <typename REALTYPE>
inline void set_aPp5_dirac(REALTYPE *yL, REALTYPE a, REALTYPE *xL,
                                                              int Nc)
{
  for(int ic = 0; ic < Nc; ++ic){
    for(int j = 0; j < VLEN; ++j){
      yL[j + VLEN*(ic + Nc* 0)] = a * ( xL[j + VLEN*(ic + Nc* 0)]
                                      + xL[j + VLEN*(ic + Nc* 2)]);
      yL[j + VLEN*(ic + Nc* 1)] = a * ( xL[j + VLEN*(ic + Nc* 1)]
                                      + xL[j + VLEN*(ic + Nc* 3)]);
      yL[j + VLEN*(ic + Nc* 2)] = a * ( xL[j + VLEN*(ic + Nc* 2)]
                                      + xL[j + VLEN*(ic + Nc* 0)]);
      yL[j + VLEN*(ic + Nc* 3)] = a * ( xL[j + VLEN*(ic + Nc* 3)]
                                      + xL[j + VLEN*(ic + Nc* 1)]);
    }
  }

}

//====================================================================
template <typename REALTYPE>
inline void set_aPm5_dirac(REALTYPE *yL, REALTYPE a, REALTYPE *xL,
                                                              int Nc)
{
  for(int ic = 0; ic < Nc; ++ic){
    for(int j = 0; j < VLEN; ++j){
      yL[j + VLEN*(ic + Nc* 0)] = a * ( xL[j + VLEN*(ic + Nc* 0)]
                                      - xL[j + VLEN*(ic + Nc* 2)]);
      yL[j + VLEN*(ic + Nc* 1)] = a * ( xL[j + VLEN*(ic + Nc* 1)]
                                      - xL[j + VLEN*(ic + Nc* 3)]);
      yL[j + VLEN*(ic + Nc* 2)] = a * ( xL[j + VLEN*(ic + Nc* 2)]
                                      - xL[j + VLEN*(ic + Nc* 0)]);
      yL[j + VLEN*(ic + Nc* 3)] = a * ( xL[j + VLEN*(ic + Nc* 3)]
                                      - xL[j + VLEN*(ic + Nc* 1)]);
    }
  }

}

//====================================================================
template <typename REALTYPE>
inline void add_aPp5_dirac(REALTYPE *yL, REALTYPE a, REALTYPE *xL,
                                                              int Nc)
{
  for(int ic = 0; ic < Nc; ++ic){
    for(int j = 0; j < VLEN; ++j){
      yL[j + VLEN*(ic + Nc* 0)] += a * ( xL[j + VLEN*(ic + Nc* 0)]
                                       + xL[j + VLEN*(ic + Nc* 2)]);
      yL[j + VLEN*(ic + Nc* 1)] += a * ( xL[j + VLEN*(ic + Nc* 1)]
                                       + xL[j + VLEN*(ic + Nc* 3)]);
      yL[j + VLEN*(ic + Nc* 2)] += a * ( xL[j + VLEN*(ic + Nc* 2)]
                                       + xL[j + VLEN*(ic + Nc* 0)]);
      yL[j + VLEN*(ic + Nc* 3)] += a * ( xL[j + VLEN*(ic + Nc* 3)]
                                       + xL[j + VLEN*(ic + Nc* 1)]);
    }
  }

}

//====================================================================
template <typename REALTYPE>
inline void add_aPm5_dirac(REALTYPE *yL, REALTYPE a, REALTYPE *xL,
                                                              int Nc)
{
  for(int ic = 0; ic < Nc; ++ic){
    for(int j = 0; j < VLEN; ++j){
      yL[j + VLEN*(ic + Nc* 0)] += a * ( xL[j + VLEN*(ic + Nc* 0)]
                                       - xL[j + VLEN*(ic + Nc* 2)]);
      yL[j + VLEN*(ic + Nc* 1)] += a * ( xL[j + VLEN*(ic + Nc* 1)]
                                       - xL[j + VLEN*(ic + Nc* 3)]);
      yL[j + VLEN*(ic + Nc* 2)] += a * ( xL[j + VLEN*(ic + Nc* 2)]
                                       - xL[j + VLEN*(ic + Nc* 0)]);
      yL[j + VLEN*(ic + Nc* 3)] += a * ( xL[j + VLEN*(ic + Nc* 3)]
                                       - xL[j + VLEN*(ic + Nc* 1)]);
    }
  }

}

//====================================================================
template <typename REALTYPE>
inline void set_aPp5_chiral(REALTYPE *yL, REALTYPE a, REALTYPE *xL,
                                                              int Nc)
{
  for(int ic = 0; ic < Nc; ++ic){
    for(int j = 0; j < VLEN; ++j){
      yL[j + VLEN*(ic + Nc* 0)] = 2.0 * a * xL[j + VLEN*(ic + Nc* 0)];
      yL[j + VLEN*(ic + Nc* 1)] = 2.0 * a * xL[j + VLEN*(ic + Nc* 1)];
      yL[j + VLEN*(ic + Nc* 2)] = 0.0;
      yL[j + VLEN*(ic + Nc* 3)] = 0.0;
    }
  }

}

//====================================================================
template <typename REALTYPE>
inline void set_aPm5_chiral(REALTYPE *yL, REALTYPE a, REALTYPE *xL,
                                                              int Nc)
{
  for(int ic = 0; ic < Nc; ++ic){
    for(int j = 0; j < VLEN; ++j){
      yL[j + VLEN*(ic + Nc* 0)] = 0.0;
      yL[j + VLEN*(ic + Nc* 1)] = 0.0;
      yL[j + VLEN*(ic + Nc* 2)] = 2.0 * a * xL[j + VLEN*(ic + Nc* 2)];
      yL[j + VLEN*(ic + Nc* 3)] = 2.0 * a * xL[j + VLEN*(ic + Nc* 3)];
    }
  }

}

} // nameless namespace end

#endif
//====================================================================
