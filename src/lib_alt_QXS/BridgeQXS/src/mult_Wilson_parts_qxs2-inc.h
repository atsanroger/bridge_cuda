/*!
      @file    mult_Wilson_parts_qxs2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_WILSON_PARTS_QXS_2_H
#define MULT_WILSON_PARTS_QXS_2_H

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
  inline void mult_wilson_xpb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_xp(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_xp(pg, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_u(pg, ut10, ut11, ut12, ut13, ut14, ut15,
             &u[VLEN * (2 * ic)]);
      mult_uv(pg, wt1r, wt1i,
              ut10, ut11, ut12, ut13, ut14, ut15,
              vt10, vt11, vt12, vt13, vt14, vt15);
      mult_uv(pg, wt2r, wt2i,
              ut10, ut11, ut12, ut13, ut14, ut15,
              vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_xp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
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
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_xm(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_xm(pg, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic]);
      mult_udv(pg, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_xm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_yp1(REALTYPE *buf, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t vt1r, vt1i, vt2r, vt2i;
      set_sp2_yp(pg, vt1r, vt1i, vt2r, vt2i, v1, ic);
      save_vec(pg, &buf[VLEN * (2 * ic)], vt1r);
      save_vec(pg, &buf[VLEN * (2 * ic + 1)], vt1i);
      save_vec(pg, &buf[VLEN * (2 * ic + NVC)], vt2r);
      save_vec(pg, &buf[VLEN * (2 * ic + 1 + NVC)], vt2i);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_yp2(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    load_vec(pg, vt10, &buf[VLEN * 0]);
    load_vec(pg, vt11, &buf[VLEN * 1]);
    load_vec(pg, vt12, &buf[VLEN * 2]);
    load_vec(pg, vt13, &buf[VLEN * 3]);
    load_vec(pg, vt14, &buf[VLEN * 4]);
    load_vec(pg, vt15, &buf[VLEN * 5]);

    load_vec(pg, vt20, &buf[VLEN * (0 + NVC)]);
    load_vec(pg, vt21, &buf[VLEN * (1 + NVC)]);
    load_vec(pg, vt22, &buf[VLEN * (2 + NVC)]);
    load_vec(pg, vt23, &buf[VLEN * (3 + NVC)]);
    load_vec(pg, vt24, &buf[VLEN * (4 + NVC)]);
    load_vec(pg, vt25, &buf[VLEN * (5 + NVC)]);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_u(pg, ut10, ut11, ut12, ut13, ut14, ut15,
             &u[VLEN * (2 * ic)]);
      mult_uv(pg, wt1r, wt1i,
              ut10, ut11, ut12, ut13, ut14, ut15,
              vt10, vt11, vt12, vt13, vt14, vt15);
      mult_uv(pg, wt2r, wt2i,
              ut10, ut11, ut12, ut13, ut14, ut15,
              vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_yp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ypb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_yp(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_yp(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_yp(pg, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_u(pg, ut10, ut11, ut12, ut13, ut14, ut15,
             &u[VLEN * (2 * ic)]);
      mult_uv(pg, wt1r, wt1i,
              ut10, ut11, ut12, ut13, ut14, ut15,
              vt10, vt11, vt12, vt13, vt14, vt15);
      mult_uv(pg, wt2r, wt2i,
              ut10, ut11, ut12, ut13, ut14, ut15,
              vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_yp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ypb(REALTYPE *v2, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_yp(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_yp(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_yp(pg, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_u(pg, ut10, ut11, ut12, ut13, ut14, ut15,
             &u[VLEN * (2 * ic)]);
      mult_uv(pg, wt1r, wt1i,
              ut10, ut11, ut12, ut13, ut14, ut15,
              vt10, vt11, vt12, vt13, vt14, vt15);
      mult_uv(pg, wt2r, wt2i,
              ut10, ut11, ut12, ut13, ut14, ut15,
              vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_yp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ym1(REALTYPE *buf, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_ym(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_ym(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_ym(pg, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic]);
      mult_udv(pg, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);

      mult_udv(pg, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      save_vec(pg, &buf[VLEN * (2 * ic)], wt1r);
      save_vec(pg, &buf[VLEN * (2 * ic + 1)], wt1i);

      save_vec(pg, &buf[VLEN * (2 * ic + NVC)], wt2r);
      save_vec(pg, &buf[VLEN * (2 * ic + 1 + NVC)], wt2i);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ym2(Vsimd_t *v2, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t wt1r, wt1i, wt2r, wt2i;
      load_vec(pg, wt1r, &buf[VLEN * (2 * ic)]);
      load_vec(pg, wt1i, &buf[VLEN * (2 * ic + 1)]);
      load_vec(pg, wt2r, &buf[VLEN * (2 * ic + NVC)]);
      load_vec(pg, wt2i, &buf[VLEN * (2 * ic + 1 + NVC)]);
      set_sp4_ym(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ymb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_ym(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_ym(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_ym(pg, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic]);
      mult_udv(pg, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_ym(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
} // nameless namespace end

#endif
