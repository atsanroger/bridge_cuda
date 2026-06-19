/*!
        @file    fopr_Wilson_impl_SU_N-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_WILSON_IMPL_INC_INCLUDED
#define FOPR_WILSON_IMPL_INC_INCLUDED

#define  NCOL    m_Nc
#define  NVC     m_Nvc
#define  ND      m_Nd

// the following macro varable us used only in inline functions
#define  NC      Nc

#define  ID1     (0)
#define  ID2     (2 * Nc)
#define  ID3     (4 * Nc)
#define  ID4     (6 * Nc)


namespace {
  std::string imple_Nc() { return "SU(N)"; }

  void check_Nc()
  {
    vout.general(CommonParameters::Vlevel(),
                 "  Gauge group implementation: general SU(N).\n");
  }


  inline double mult_uv_r(const double *g, const double *w, const int Nc)
  {
    double a = 0.0;
    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i] * w[2 * i] - g[2 * i + 1] * w[2 * i + 1];
    }
    return a;
  }


  inline double mult_uv_i(const double *g, const double *w, const int Nc)
  {
    double a = 0.0;
    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i] * w[2 * i + 1] + g[2 * i + 1] * w[2 * i];
    }
    return a;
  }


  inline double mult_udagv_r(const double *g, const double *w, const int Nc)
  {
    double a = 0.0;
    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i * Nc] * w[2 * i] + g[2 * i * Nc + 1] * w[2 * i + 1];
    }
    return a;
  }


  inline double mult_udagv_i(const double *g, const double *w, const int Nc)
  {
    double a = 0.0;
    for (int i = 0; i < Nc; ++i) {
      a += g[2 * i * Nc] * w[2 * i + 1] - g[2 * i * Nc + 1] * w[2 * i];
    }
    return a;
  }
} // end of nameless namespace

#endif
//============================================================END=====
