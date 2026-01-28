/*!
        @file    fopr_Wilson_impl_common-inc.h
        @brief
        @author  HIDeo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_WILSON_IMPL_COMMON_INC_INCLUDED
#define FOPR_WILSON_IMPL_COMMON_INC_INCLUDED


namespace {
  inline void mult_gamma5_dirac(double *v, const double *w, const int Nc)
  {
    for (int ivc = 0; ivc < 2 * NC; ++ivc) {
      v[ivc + ID1] = w[ivc + ID3];
      v[ivc + ID2] = w[ivc + ID4];
      v[ivc + ID3] = w[ivc + ID1];
      v[ivc + ID4] = w[ivc + ID2];
    }
  }


  inline void mult_gamma5_chiral(double *v, const double *w, const int Nc)
  {
    for (int ivc = 0; ivc < 2 * NC; ++ivc) {
      v[ivc + ID1] = w[ivc + ID1];
      v[ivc + ID2] = w[ivc + ID2];
      v[ivc + ID3] = -w[ivc + ID3];
      v[ivc + ID4] = -w[ivc + ID4];
    }
  }


  inline void set_sp2_xp(double *vt1, double *vt2,
                         const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = w[icr + ID1] - w[ici + ID4];
      vt1[ici] = w[ici + ID1] + w[icr + ID4];
      vt2[icr] = w[icr + ID2] - w[ici + ID3];
      vt2[ici] = w[ici + ID2] + w[icr + ID3];
    }
  }


  inline void set_sp4_xp(double *v, double w1r, double w1i,
                         double w2r, double w2i, const int Nc)
  {
    v[0 + ID1] += w1r;
    v[1 + ID1] += w1i;
    v[0 + ID2] += w2r;
    v[1 + ID2] += w2i;
    v[0 + ID3] += w2i;
    v[1 + ID3] -= w2r;
    v[0 + ID4] += w1i;
    v[1 + ID4] -= w1r;
  }


  inline void set_sp2_xm(double *vt1, double *vt2,
                         const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = w[icr + ID1] + w[ici + ID4];
      vt1[ici] = w[ici + ID1] - w[icr + ID4];
      vt2[icr] = w[icr + ID2] + w[ici + ID3];
      vt2[ici] = w[ici + ID2] - w[icr + ID3];
    }
  }


  inline void set_sp4_xm(double *v, double w1r, double w1i,
                         double w2r, double w2i, const int Nc)
  {
    v[0 + ID1] += w1r;
    v[1 + ID1] += w1i;
    v[0 + ID2] += w2r;
    v[1 + ID2] += w2i;
    v[0 + ID3] -= w2i;
    v[1 + ID3] += w2r;
    v[0 + ID4] -= w1i;
    v[1 + ID4] += w1r;
  }


  inline void set_sp2_yp(double *vt1, double *vt2,
                         const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = w[icr + ID1] + w[icr + ID4];
      vt1[ici] = w[ici + ID1] + w[ici + ID4];
      vt2[icr] = w[icr + ID2] - w[icr + ID3];
      vt2[ici] = w[ici + ID2] - w[ici + ID3];
    }
  }


  inline void set_sp4_yp(double *v, double w1r, double w1i,
                         double w2r, double w2i, const int Nc)
  {
    v[0 + ID1] += w1r;
    v[1 + ID1] += w1i;
    v[0 + ID2] += w2r;
    v[1 + ID2] += w2i;
    v[0 + ID3] -= w2r;
    v[1 + ID3] -= w2i;
    v[0 + ID4] += w1r;
    v[1 + ID4] += w1i;
  }


  inline void set_sp2_ym(double *vt1, double *vt2,
                         const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = w[icr + ID1] - w[icr + ID4];
      vt1[ici] = w[ici + ID1] - w[ici + ID4];
      vt2[icr] = w[icr + ID2] + w[icr + ID3];
      vt2[ici] = w[ici + ID2] + w[ici + ID3];
    }
  }


  inline void set_sp4_ym(double *v, double w1r, double w1i,
                         double w2r, double w2i, const int Nc)
  {
    v[0 + ID1] += w1r;
    v[1 + ID1] += w1i;
    v[0 + ID2] += w2r;
    v[1 + ID2] += w2i;
    v[0 + ID3] += w2r;
    v[1 + ID3] += w2i;
    v[0 + ID4] -= w1r;
    v[1 + ID4] -= w1i;
  }


  inline void set_sp2_zp(double *vt1, double *vt2,
                         const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = w[icr + ID1] - w[ici + ID3];
      vt1[ici] = w[ici + ID1] + w[icr + ID3];
      vt2[icr] = w[icr + ID2] + w[ici + ID4];
      vt2[ici] = w[ici + ID2] - w[icr + ID4];
    }
  }


  inline void set_sp4_zp(double *v, double w1r, double w1i,
                         double w2r, double w2i, const int Nc)
  {
    v[0 + ID1] += w1r;
    v[1 + ID1] += w1i;
    v[0 + ID2] += w2r;
    v[1 + ID2] += w2i;
    v[0 + ID3] += w1i;
    v[1 + ID3] -= w1r;
    v[0 + ID4] -= w2i;
    v[1 + ID4] += w2r;
  }


  inline void set_sp2_zm(double *vt1, double *vt2,
                         const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = w[icr + ID1] + w[ici + ID3];
      vt1[ici] = w[ici + ID1] - w[icr + ID3];
      vt2[icr] = w[icr + ID2] - w[ici + ID4];
      vt2[ici] = w[ici + ID2] + w[icr + ID4];
    }
  }


  inline void set_sp4_zm(double *v, double w1r, double w1i,
                         double w2r, double w2i, const int Nc)
  {
    v[0 + ID1] += w1r;
    v[1 + ID1] += w1i;
    v[0 + ID2] += w2r;
    v[1 + ID2] += w2i;
    v[0 + ID3] -= w1i;
    v[1 + ID3] += w1r;
    v[0 + ID4] += w2i;
    v[1 + ID4] -= w2r;
  }


  inline void set_sp2_tp_dirac(double *vt1, double *vt2,
                               const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = 2.0 * w[icr + ID3];
      vt1[ici] = 2.0 * w[ici + ID3];
      vt2[icr] = 2.0 * w[icr + ID4];
      vt2[ici] = 2.0 * w[ici + ID4];
    }
  }


  inline void set_sp4_tp_dirac(double *v,
                               double w1r, double w1i,
                               double w2r, double w2i, const int Nc)
  {
    v[0 + ID3] += w1r;
    v[1 + ID3] += w1i;
    v[0 + ID4] += w2r;
    v[1 + ID4] += w2i;
  }


  inline void set_sp2_tm_dirac(double *vt1, double *vt2,
                               const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = 2.0 * w[icr + ID1];
      vt1[ici] = 2.0 * w[ici + ID1];
      vt2[icr] = 2.0 * w[icr + ID2];
      vt2[ici] = 2.0 * w[ici + ID2];
    }
  }


  inline void set_sp4_tm_dirac(double *v,
                               double w1r, double w1i,
                               double w2r, double w2i, const int Nc)
  {
    v[0 + ID1] += w1r;
    v[1 + ID1] += w1i;
    v[0 + ID2] += w2r;
    v[1 + ID2] += w2i;
  }


  inline void set_sp2_tp_chiral(double *vt1, double *vt2,
                                const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = w[icr + ID1] + w[icr + ID3];
      vt1[ici] = w[ici + ID1] + w[ici + ID3];
      vt2[icr] = w[icr + ID2] + w[icr + ID4];
      vt2[ici] = w[ici + ID2] + w[ici + ID4];
    }
  }


  inline void set_sp4_tp_chiral(double *v,
                                double w1r, double w1i,
                                double w2r, double w2i, const int Nc)
  {
    v[0 + ID1] += w1r;
    v[1 + ID1] += w1i;
    v[0 + ID2] += w2r;
    v[1 + ID2] += w2i;
    v[0 + ID3] += w1r;
    v[1 + ID3] += w1i;
    v[0 + ID4] += w2r;
    v[1 + ID4] += w2i;
  }


  inline void set_sp2_tm_chiral(double *vt1, double *vt2,
                                const double *w, const int Nc)
  {
    for (int ic = 0; ic < NC; ++ic) {
      int icr = 2 * ic;
      int ici = 2 * ic + 1;
      vt1[icr] = w[icr + ID1] - w[icr + ID3];
      vt1[ici] = w[ici + ID1] - w[ici + ID3];
      vt2[icr] = w[icr + ID2] - w[icr + ID4];
      vt2[ici] = w[ici + ID2] - w[ici + ID4];
    }
  }


  inline void set_sp4_tm_chiral(double *v,
                                double w1r, double w1i,
                                double w2r, double w2i, const int Nc)
  {
    v[0 + ID1] += w1r;
    v[1 + ID1] += w1i;
    v[0 + ID2] += w2r;
    v[1 + ID2] += w2i;
    v[0 + ID3] -= w1r;
    v[1 + ID3] -= w1i;
    v[0 + ID4] -= w2r;
    v[1 + ID4] -= w2i;
  }
} // end of nameless namespace

#endif
//============================================================END=====
