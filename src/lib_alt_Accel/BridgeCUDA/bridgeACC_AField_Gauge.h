/*!
      @file    bridgeACC_AField_Gauge.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
      @version $LastChangedRevision: 2535 $
*/

#ifndef BRIDGEACC_AFIELD_GAUGE_INCLUDED
#define BRIDGEACC_AFIELD_GAUGE_INCLUDED

namespace BridgeACC {

// real_t = double

void multadd_Gnn(double* u, const int exu,
                 double* v, const int exv,
                 double* w, const int exw,
                 const double a, const int nst);

void mult_Gnn(double* u, const int exu,
              double* v, const int exv,
              double* w, const int exw, const int nst);

void multadd_Gnd(double* u, const int exu,
                 double* v, const int exv,
                 double* w, const int exw,
                 const double a, const int nst);

void mult_Gnd(double* u, const int exu,
              double* v, const int exv,
              double* w, const int exw, const int nst);

void multadd_Gdn(double* u, const int exu,
                 double* v, const int exv,
                 double* w, const int exw,
		 const double a, const int nst);

void mult_Gdn(double* u, const int exu,
              double* v, const int exv,
              double* w, const int exw, const int nst);

void mult_Gdd(double* u, const int exu,
              double* v, const int exv,
              double* w, const int exw, const int nst);

void ah_G(double *u, const int ex, const int nst);

void at_G(double *u, const int ex, const int nst);

void add_unit(double *u, const int ex, const double a, const int nst);

void inverse_dag(double* uinv, const int ex1,
                 double* u, const int ex2,
                 const int nst);

// real_t = float

void multadd_Gnn(float* u, const int exu,
                 float* v, const int exv,
                 float* w, const int exw,
                 const float a, const int nst);

void mult_Gnn(float* u, const int exu,
              float* v, const int exv,
              float* w, const int exw, const int nst);

void multadd_Gnd(float* u, const int exu,
                 float* v, const int exv,
                 float* w, const int exw,
                 const float a, const int nst);

void mult_Gnd(float* u, const int exu,
              float* v, const int exv,
              float* w, const int exw, const int nst);

void multadd_Gdn(float* u, const int exu,
                 float* v, const int exv,
                 float* w, const int exw,
		 const float a, const int nst);

void mult_Gdn(float* u, const int exu,
              float* v, const int exv,
              float* w, const int exw, const int nst);

void mult_Gdd(float* u, const int exu,
              float* v, const int exv,
              float* w, const int exw, const int nst);

void ah_G(float *u, const int ex, const int nst);

void at_G(float *u, const int ex, const int nst);

void add_unit(float *u, const int ex, const float a, const int nst);

void inverse_dag(float* uinv, const int ex1,
                 float* u, const int ex2,
                 const int nst);

}

//============================================================END=====
#endif
