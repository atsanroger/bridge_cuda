/*!
        @file    fopr_Wilson_impl_SU2-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_WILSON_IMPL_INC_INCLUDED
#define FOPR_WILSON_IMPL_INC_INCLUDED

#define  NC      2
#define  NCOL    2
#define  NVC     4
#define  ND      4

#define  ID1     0
#define  ID2     4
#define  ID3     8
#define  ID4     12

namespace {
  std::string imple_Nc() { return "SU(2)"; }

  void check_Nc()
  {
    vout.general(CommonParameters::Vlevel(),
                 "  Gauge group implementation: SU(2).\n");
  }


  inline double mult_uv_r(const double *u, const double *v, const int Nc)
  {
    return u[0] * v[0] - u[1] * v[1]
           + u[2] * v[2] - u[3] * v[3];
  }


  inline double mult_uv_i(const double *u, const double *v, const int Nc)
  {
    return u[0] * v[1] + u[1] * v[0]
           + u[2] * v[3] + u[3] * v[2];
  }


  inline double mult_udagv_r(const double *u, const double *v, const int Nc)
  {
    return u[0] * v[0] + u[1] * v[1]
           + u[4] * v[2] + u[5] * v[3];
  }


  inline double mult_udagv_i(const double *u, const double *v, const int Nc)
  {
    return u[0] * v[1] - u[1] * v[0]
           + u[4] * v[3] - u[5] * v[2];
  }
} // end of nameless namespace

#endif
//============================================================END=====
