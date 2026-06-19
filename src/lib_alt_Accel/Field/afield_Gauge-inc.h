/*!
      @file    afield_Gauge-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-09-10 17:39:07 #$
      @version $LastChangedRevision: 2537 $
*/

#ifndef ACCEL_AFIELD_GAUGE_INC_INCLUDED
#define ACCEL_AFIELD_GAUGE_INC_INCLUDED

//#include <cstdlib> 

#include "lib_alt_Accel/BridgeACC/bridgeACC_AField_Gauge.h"

//====================================================================
template <typename REALTYPE>
void multadd_Gnn(AField<REALTYPE,ACCEL>& u, const int exu,
                 const AField<REALTYPE,ACCEL>& v, const int exv,
                 const AField<REALTYPE,ACCEL>& w, const int exw,
                 const REALTYPE a)
{
  typedef AField<REALTYPE,ACCEL> AFIELD;
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = v.nvol_pad();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(ith == 0)
    BridgeACC::multadd_Gnn(up, exu, vp, exv, wp, exw, a, Nst);

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void mult_Gnn(AField<REALTYPE,ACCEL>& u, const int exu,
              const AField<REALTYPE,ACCEL>& v, const int exv,
              const AField<REALTYPE,ACCEL>& w, const int exw)
{
  typedef AField<REALTYPE,ACCEL> AFIELD;
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = v.nvol_pad();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(ith == 0) BridgeACC::mult_Gnn(up, exu, vp, exv, wp, exw, Nst);

#pragma omp barrier

}

//====================================================================
template <typename REALTYPE>
void multadd_Gnd(AField<REALTYPE,ACCEL>& u, const int exu,
                 const AField<REALTYPE,ACCEL>& v, const int exv,
                 const AField<REALTYPE,ACCEL>& w, const int exw,
                 const REALTYPE a)
{
  typedef AField<REALTYPE,ACCEL> AFIELD;
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = v.nvol_pad();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(ith == 0)
    BridgeACC::multadd_Gnd(up, exu, vp, exv, wp, exw, a, Nst);

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void mult_Gnd(AField<REALTYPE,ACCEL>& u, const int exu,
              const AField<REALTYPE,ACCEL>& v, const int exv,
              const AField<REALTYPE,ACCEL>& w, const int exw)
{
  typedef AField<REALTYPE,ACCEL> AFIELD;
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = v.nvol_pad();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(ith == 0) BridgeACC::mult_Gnd(up, exu, vp, exv, wp, exw, Nst);

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void multadd_Gdn(AField<REALTYPE,ACCEL>& u, const int exu,
                 const AField<REALTYPE,ACCEL>& v, const int exv,
                 const AField<REALTYPE,ACCEL>& w, const int exw,
                 const REALTYPE a)
{
  typedef AField<REALTYPE,ACCEL> AFIELD;
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = v.nvol_pad();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(ith == 0)
    BridgeACC::multadd_Gdn(up, exu, vp, exv, wp, exw, a, Nst);

#pragma omp barrier

}

//====================================================================
template <typename REALTYPE>
void mult_Gdn(AField<REALTYPE,ACCEL>& u, const int exu,
              const AField<REALTYPE,ACCEL>& v, const int exv,
              const AField<REALTYPE,ACCEL>& w, const int exw)
{
  typedef AField<REALTYPE,ACCEL> AFIELD;
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = v.nvol_pad();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(ith == 0) BridgeACC::mult_Gdn(up, exu, vp, exv, wp, exw, Nst);

#pragma omp barrier

}

//====================================================================
template <typename REALTYPE>
void mult_Gdd(AField<REALTYPE,ACCEL>& u, const int exu,
              const AField<REALTYPE,ACCEL>& v, const int exv,
              const AField<REALTYPE,ACCEL>& w, const int exw)
{
  typedef AField<REALTYPE,ACCEL> AFIELD;
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = v.nvol_pad();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  if(ith == 0) BridgeACC::mult_Gdd(up, exu, vp, exv, wp, exw, Nst);

#pragma omp barrier
}

//====================================================================
// anti-hermitian
template <typename REALTYPE>
void ah_G(AField<REALTYPE,ACCEL>& u, const int ex)
{
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = u.nvol();
  real_t* up = u.ptr(0);

  if(ith == 0) BridgeACC::ah_G(up, ex, Nst);

#pragma omp barrier

}

//====================================================================
// anti-hermitian traceless
template <typename REALTYPE>
void at_G(AField<REALTYPE,ACCEL>& u, const int ex)
{
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = u.nvol();
  real_t* up = u.ptr(0);

  if(ith == 0) BridgeACC::at_G(up, ex, Nst);

#pragma omp barrier

}

//====================================================================
template <typename REALTYPE>
void add_unit(AField<REALTYPE,ACCEL>& u, const int ex, REALTYPE a)
{         // u = u + a * I (I: unit matrix)

  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  int Nst = u.nvol();
  real_t* up = u.ptr(0);

  if(ith == 0) BridgeACC::add_unit(up, ex, a, Nst);

#pragma omp barrier

}

//============================================================END=====
#endif
