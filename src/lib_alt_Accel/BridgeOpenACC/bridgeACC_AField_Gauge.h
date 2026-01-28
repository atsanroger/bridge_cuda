/*!
      @file    bridgeACC_AField_Gauge.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef BRIDGEACC_AFIELD_GAUGE_INCLUDED
#define BRIDGEACC_AFIELD_GAUGE_INCLUDED

namespace BridgeACC {

// real_t = double

void multadd_Gnn(double *restrict u, const int exu,
                 double *restrict v, const int exv,
                 double *restrict w, const int exw,
                 const double a, const int nst);

void mult_Gnn(double *restrict u, const int exu,
              double *restrict v, const int exv,
              double *restrict w, const int exw, const int nst);

void multadd_Gnd(double *restrict u, const int exu,
                 double *restrict v, const int exv,
                 double *restrict w, const int exw,
                 const double a, const int nst);

void mult_Gnd(double *restrict u, const int exu,
              double *restrict v, const int exv,
              double *restrict w, const int exw, const int nst);

void multadd_Gdn(double *restrict u, const int exu,
                 double *restrict v, const int exv,
                 double *restrict w, const int exw,
		 const double a, const int nst);

void mult_Gdn(double *restrict u, const int exu,
              double *restrict v, const int exv,
              double *restrict w, const int exw, const int nst);

void mult_Gdd(double *restrict u, const int exu,
              double *restrict v, const int exv,
              double *restrict w, const int exw, const int nst);

void ah_G(double *u, const int ex, const int nst);

void at_G(double *u, const int ex, const int nst);

void add_unit(double *u, const int ex, const double a, const int nst);

void inverse_dag(double *restrict uinv, const int ex1,
                 const double *restrict u, const int ex2,
                 const int nst);

// real_t = float

void multadd_Gnn(float *restrict u, const int exu,
                 float *restrict v, const int exv,
                 float *restrict w, const int exw,
                 const float a, const int nst);

void mult_Gnn(float *restrict u, const int exu,
              float *restrict v, const int exv,
              float *restrict w, const int exw, const int nst);

void multadd_Gnd(float *restrict u, const int exu,
                 float *restrict v, const int exv,
                 float *restrict w, const int exw,
                 const float a, const int nst);

void mult_Gnd(float *restrict u, const int exu,
              float *restrict v, const int exv,
              float *restrict w, const int exw, const int nst);

void multadd_Gdn(float *restrict u, const int exu,
                 float *restrict v, const int exv,
                 float *restrict w, const int exw,
		 const float a, const int nst);

void mult_Gdn(float *restrict u, const int exu,
              float *restrict v, const int exv,
              float *restrict w, const int exw, const int nst);

void mult_Gdd(float *restrict u, const int exu,
              float *restrict v, const int exv,
              float *restrict w, const int exw, const int nst);

void ah_G(float *u, const int ex, const int nst);

void at_G(float *u, const int ex, const int nst);

void add_unit(float *u, const int ex, const float a, const int nst);

void inverse_dag(float *restrict uinv, const int ex1,
                 const float *restrict u, const int ex2,
                 const int nst);

}

//============================================================END=====
#endif
