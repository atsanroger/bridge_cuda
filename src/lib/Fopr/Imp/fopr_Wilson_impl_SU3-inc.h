/*!
        @file    fopr_Wilson_impl_SU3-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_WILSON_IMPL_INC_INCLUDED
#define FOPR_WILSON_IMPL_INC_INCLUDED

#define  NCOL    3
#define  NVC     6
#define  ND      4

// the following macro varable us used only in inline functions
#define  NC      3

#define  ID1     0
#define  ID2     6
#define  ID3     12
#define  ID4     18

namespace {
  std::string imple_Nc() { return "SU(3)"; }

  void check_Nc()
  {
    vout.general(CommonParameters::Vlevel(),
                 "  Gauge group implementation: SU(3).\n");
  }


  inline double mult_uv_r(const double *g, const double *w, const int Nc)
  {
    return g[0] * w[0] - g[1] * w[1]
           + g[2] * w[2] - g[3] * w[3]
           + g[4] * w[4] - g[5] * w[5];
  }


  inline double mult_uv_i(const double *g, const double *w, const int Nc)
  {
    return g[0] * w[1] + g[1] * w[0]
           + g[2] * w[3] + g[3] * w[2]
           + g[4] * w[5] + g[5] * w[4];
  }


  inline double mult_udagv_r(const double *g, const double *w, const int Nc)
  {
    return g[0] * w[0] + g[1] * w[1]
           + g[6] * w[2] + g[7] * w[3]
           + g[12] * w[4] + g[13] * w[5];
  }


  inline double mult_udagv_i(const double *g, const double *w, const int Nc)
  {
    return g[0] * w[1] - g[1] * w[0]
           + g[6] * w[3] - g[7] * w[2]
           + g[12] * w[5] - g[13] * w[4];
  }
} // end of nameless namespace

#endif
//============================================================END=====
