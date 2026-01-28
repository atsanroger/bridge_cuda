/*!
      @file    mult_Wilson_simd2_parts-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
       @date   $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_WILSON_QXS_PARTS_H
#define MULT_WILSON_QXS_PARTS_H

namespace {
//====================================================================
  inline void check_setup()
  {
    /*
    if(VLEN2 < 2){
      vout.crucial("VLEN2 = %d is too small for this implementation\n",
                   VLEN2);
      exit(EXIT_FAILURE);
    }
    */
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xp1(REALTYPE *buf, REALTYPE *v1)
  {
    REALTYPE vt[NVCD];

    load_vec1(vt, v1, 0, NVCD);
    set_sp2_xp1(buf, vt);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xp2(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    shift_vec1_bw(vt1, &buf[0], NVC);
    shift_vec1_bw(vt2, &buf[NVC], NVC);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      mult_uv(wt1, &ut[2 * ic], vt1, NC);
      mult_uv(wt2, &ut[2 * ic], vt2, NC);
      set_sp4_xp(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xpb(Vsimd_t *v2,
                              REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_xp(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      mult_uv(wt1, &ut[2 * ic], vt1, NC);
      mult_uv(wt2, &ut[2 * ic], vt2, NC);
      set_sp4_xp(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xm1(REALTYPE *buf, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[VLEN * NVC], vt2[VLEN * NVC];
    set_sp2_xm(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[NVC], wt2[NVC];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = NVC * ic;
      mult_udagv(&wt1[2 * ic], &ut[ic2], vt1, NC);
      mult_udagv(&wt2[2 * ic], &ut[ic2], vt2, NC);
    }

    for (int ic = 0; ic < NC; ++ic) {
      save_vec1(&buf[0], wt1, VLEN - 1, NVC);
      save_vec1(&buf[NVC], wt2, VLEN - 1, NVC);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xm2(Vsimd_t *v2, REALTYPE *buf)
  {
    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      shift_vec1_fw(wt1, &buf[2 * ic], 2);
      shift_vec1_fw(wt2, &buf[2 * ic + NVC], 2);
      set_sp4_xm(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xmb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_xm(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = NVC * ic;
      int ic3 = ND * 2 * ic;
      mult_udagv(wt1, &ut[ic2], vt1, NC);
      mult_udagv(wt2, &ut[ic2], vt2, NC);
      set_sp4_xm(&v2[ic3], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_yp1(REALTYPE *buf, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_yp(vt1, vt2, v1);

    save_vec(&buf[0], vt1, NVC);
    save_vec(&buf[VLEN * NVC], vt2, NVC);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_yp2(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    load_vec(vt1, &buf[0], NVC);
    load_vec(vt2, &buf[VLEN * NVC], NVC);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      mult_uv(wt1, &ut[2 * ic], vt1, NC);
      mult_uv(wt2, &ut[2 * ic], vt2, NC);
      set_sp4_yp(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ypb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_yp(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      mult_uv(wt1, &ut[2 * ic], vt1, NC);
      mult_uv(wt2, &ut[2 * ic], vt2, NC);
      set_sp4_yp(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ym1(REALTYPE *buf, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_ym(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[NVC], wt2[NVC];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = NVC * ic;
      mult_udagv(&wt1[2 * ic], &ut[ic2], vt1, NC);
      mult_udagv(&wt2[2 * ic], &ut[ic2], vt2, NC);
    }

    save_vec(&buf[0], wt1, NVC);
    save_vec(&buf[VLEN * NVC], wt2, NVC);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ym2(Vsimd_t *v2, REALTYPE *buf)
  {
    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      load_vec(wt1, &buf[VLEN * 2 * ic], 2);
      load_vec(wt2, &buf[VLEN * 2 * (ic + NC)], 2);
      set_sp4_ym(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ymb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_ym(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = NVC * ic;
      int ic3 = ND * 2 * ic;
      mult_udagv(wt1, &ut[ic2], vt1, NC);
      mult_udagv(wt2, &ut[ic2], vt2, NC);
      set_sp4_ym(&v2[ic3], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zp1(REALTYPE *buf, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_zp(vt1, vt2, v1);

    save_vec(&buf[0], vt1, NVC);
    save_vec(&buf[VLEN * NVC], vt2, NVC);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zp2(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    load_vec(vt1, &buf[0], NVC);
    load_vec(vt2, &buf[VLEN * NVC], NVC);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      mult_uv(wt1, &ut[2 * ic], vt1, NC);
      mult_uv(wt2, &ut[2 * ic], vt2, NC);
      set_sp4_zp(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zpb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_zp(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      mult_uv(wt1, &ut[2 * ic], vt1, NC);
      mult_uv(wt2, &ut[2 * ic], vt2, NC);
      set_sp4_zp(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zm1(REALTYPE *buf, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_zm(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[NVC], wt2[NVC];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = NVC * ic;
      mult_udagv(&wt1[2 * ic], &ut[ic2], vt1, NC);
      mult_udagv(&wt2[2 * ic], &ut[ic2], vt2, NC);
    }

    save_vec(&buf[0], wt1, NVC);
    save_vec(&buf[VLEN * NVC], wt2, NVC);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zm2(Vsimd_t *v2, REALTYPE *buf)
  {
    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      load_vec(wt1, &buf[VLEN * 2 * ic], 2);
      load_vec(wt2, &buf[VLEN * 2 * (ic + NC)], 2);
      set_sp4_zm(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zmb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_zm(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = NVC * ic;
      int ic3 = ND * 2 * ic;
      mult_udagv(wt1, &ut[ic2], vt1, NC);
      mult_udagv(wt2, &ut[ic2], vt2, NC);
      set_sp4_zm(&v2[ic3], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tp1_dirac(REALTYPE *buf, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_tp_dirac(vt1, vt2, v1);

    save_vec(&buf[0], vt1, NVC);
    save_vec(&buf[VLEN * NVC], vt2, NVC);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tp2_dirac(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    load_vec(vt1, &buf[0], NVC);
    load_vec(vt2, &buf[VLEN * NVC], NVC);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      mult_uv(wt1, &ut[2 * ic], vt1, NC);
      mult_uv(wt2, &ut[2 * ic], vt2, NC);
      set_sp4_tp_dirac(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tpb_dirac(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_tp_dirac(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      mult_uv(wt1, &ut[2 * ic], vt1, NC);
      mult_uv(wt2, &ut[2 * ic], vt2, NC);
      set_sp4_tp_dirac(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tm1_dirac(REALTYPE *buf, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_tm_dirac(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[NVC], wt2[NVC];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = NVC * ic;
      mult_udagv(&wt1[2 * ic], &ut[ic2], vt1, NC);
      mult_udagv(&wt2[2 * ic], &ut[ic2], vt2, NC);
    }

    save_vec(&buf[0], wt1, NVC);
    save_vec(&buf[VLEN * NVC], wt2, NVC);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tm2_dirac(Vsimd_t *v2, REALTYPE *buf)
  {
    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = ND * 2 * ic;
      load_vec(wt1, &buf[VLEN * 2 * ic], 2);
      load_vec(wt2, &buf[VLEN * 2 * (ic + NC)], 2);
      set_sp4_tm_dirac(&v2[ic2], wt1, wt2);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tmb_dirac(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    Vsimd_t vt1[NVC], vt2[NVC];
    set_sp2_tm_dirac(vt1, vt2, v1);

    Vsimd_t ut[NDF];
    load_vec(ut, u, NDF);

    Vsimd_t wt1[2], wt2[2];
    for (int ic = 0; ic < NC; ++ic) {
      int ic2 = NVC * ic;
      int ic3 = ND * 2 * ic;
      mult_udagv(wt1, &ut[ic2], vt1, NC);
      mult_udagv(wt2, &ut[ic2], vt2, NC);
      set_sp4_tm_dirac(&v2[ic3], wt1, wt2);
    }
  }


//====================================================================
} // nameless namespace end

#endif
