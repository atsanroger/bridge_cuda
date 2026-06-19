/*!
      @file    mult_Wilson_eo_parts_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_WILSON_EO_PARTS_QXS_H
#define MULT_WILSON_EO_PARTS_QXS_H

namespace {
//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xp1(svbool_t& pg2, svint_t& svidx,
                                 REALTYPE *__restrict buf,
                                 REALTYPE *__restrict v1)
  {
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
    svreal_t v1r, v1i, v2r, v2i;

    for (int ic = 0; ic < NC; ++ic) {
      int icr = ND * 2 * ic;
      int ici = ND * 2 * ic + 1;

      load_vec(pg2, w1r, &v1[VLEN * (icr + ID1)]);
      load_vec(pg2, w1i, &v1[VLEN * (ici + ID1)]);
      load_vec(pg2, w2r, &v1[VLEN * (icr + ID2)]);
      load_vec(pg2, w2i, &v1[VLEN * (ici + ID2)]);
      load_vec(pg2, w3r, &v1[VLEN * (icr + ID3)]);
      load_vec(pg2, w3i, &v1[VLEN * (ici + ID3)]);
      load_vec(pg2, w4r, &v1[VLEN * (icr + ID4)]);
      load_vec(pg2, w4i, &v1[VLEN * (ici + ID4)]);

      add_vec(pg2, v1r, w1r, w4i);
      sub_vec(pg2, v1i, w1i, w4r);
      add_vec(pg2, v2r, w2r, w3i);
      sub_vec(pg2, v2i, w2i, w3r);

      int skip = (VLENY + 1) / 2;
      save_vec_scatter(pg2, &buf[skip * (2 * ic)], v1r, svidx);
      save_vec_scatter(pg2, &buf[skip * (2 * ic + 1)], v1i, svidx);
      save_vec_scatter(pg2, &buf[skip * (2 * ic + NVC)], v2r, svidx);
      save_vec_scatter(pg2, &buf[skip * (2 * ic + 1 + NVC)], v2i, svidx);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xp1(svbool_t& pg2,
                                 REALTYPE *__restrict buf,
                                 REALTYPE *__restrict v1)
  {
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
    svreal_t v1r, v1i, v2r, v2i;

    for (int ic = 0; ic < NC; ++ic) {
      int icr = ND * 2 * ic;
      int ici = ND * 2 * ic + 1;

      load_vec(pg2, w1r, &v1[VLEN * (icr + ID1)]);
      load_vec(pg2, w1i, &v1[VLEN * (ici + ID1)]);
      load_vec(pg2, w2r, &v1[VLEN * (icr + ID2)]);
      load_vec(pg2, w2i, &v1[VLEN * (ici + ID2)]);
      load_vec(pg2, w3r, &v1[VLEN * (icr + ID3)]);
      load_vec(pg2, w3i, &v1[VLEN * (ici + ID3)]);
      load_vec(pg2, w4r, &v1[VLEN * (icr + ID4)]);
      load_vec(pg2, w4i, &v1[VLEN * (ici + ID4)]);

      add_vec(pg2, v1r, w1r, w4i);
      sub_vec(pg2, v1i, w1i, w4r);
      add_vec(pg2, v2r, w2r, w3i);
      sub_vec(pg2, v2i, w2i, w3r);

      v1r = compact_vec(pg2, v1r);
      v1i = compact_vec(pg2, v1i);
      v2r = compact_vec(pg2, v2r);
      v2i = compact_vec(pg2, v2i);

      int      skip = (VLENY + 1) / 2;
      svbool_t pg1  = set_predicate_whilelt(skip);
      save_vec(pg1, &buf[skip * (2 * ic)], v1r);
      save_vec(pg1, &buf[skip * (2 * ic + 1)], v1i);
      save_vec(pg1, &buf[skip * (2 * ic + NVC)], v2r);
      save_vec(pg1, &buf[skip * (2 * ic + 1 + NVC)], v2i);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_xp2(svbool_t& pg, svbool_t& pg1, svbool_t& pg2,
                          svbool_t& pg3,
                          svreal_t& vt1r, svreal_t& vt1i,
                          svreal_t& vt2r, svreal_t& vt2i,
                          REALTYPE *__restrict v,
                          REALTYPE *__restrict buf,
                          svint_t& index, int ic)
  {
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;

    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
    load_vec(pg3, w1r, &v[VLEN * (icr + ID1)]);
    load_vec(pg3, w1i, &v[VLEN * (ici + ID1)]);
    load_vec(pg3, w2r, &v[VLEN * (icr + ID2)]);
    load_vec(pg3, w2i, &v[VLEN * (ici + ID2)]);

    load_add(pg1, w1r, &v[VLEN * (icr + ID1) + 1]);
    load_add(pg1, w1i, &v[VLEN * (ici + ID1) + 1]);
    load_add(pg1, w2r, &v[VLEN * (icr + ID2) + 1]);
    load_add(pg1, w2i, &v[VLEN * (ici + ID2) + 1]);

    load_vec(pg3, w3r, &v[VLEN * (icr + ID3)]);
    load_vec(pg3, w3i, &v[VLEN * (ici + ID3)]);
    load_vec(pg3, w4r, &v[VLEN * (icr + ID4)]);
    load_vec(pg3, w4i, &v[VLEN * (ici + ID4)]);

    svbool_t pg13 = sveor_z(pg, pg1, pg3);

    load_add(pg1, w3r, &v[VLEN * (icr + ID3) + 1]);
    load_add(pg1, w3i, &v[VLEN * (ici + ID3) + 1]);
    load_add(pg1, w4r, &v[VLEN * (icr + ID4) + 1]);
    load_add(pg1, w4i, &v[VLEN * (ici + ID4) + 1]);

    add_vec(pg13, vt1r, w1r, w4i);
    sub_vec(pg13, vt1i, w1i, w4r);
    add_vec(pg13, vt2r, w2r, w3i);
    sub_vec(pg13, vt2i, w2i, w3r);

    int skip = (VLENY + 1) / 2;
    load_add_gather(pg2, vt1r, &buf[skip * (2 * ic)], index);
    load_add_gather(pg2, vt1i, &buf[skip * (2 * ic + 1)], index);
    load_add_gather(pg2, vt2r, &buf[skip * (2 * ic + NVC)], index);
    load_add_gather(pg2, vt2i, &buf[skip * (2 * ic + 1 + NVC)], index);
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_xp2(svbool_t& pg, svbool_t& pg1, svbool_t& pg2,
                          svbool_t& pg3,
                          svreal_t& vt1r, svreal_t& vt1i,
                          svreal_t& vt2r, svreal_t& vt2i,
                          REALTYPE *__restrict v,
                          REALTYPE *__restrict buf,
                          svuint_t& idx, int ic)
  {
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;

    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
    load_vec(pg3, w1r, &v[VLEN * (icr + ID1)]);
    load_vec(pg3, w1i, &v[VLEN * (ici + ID1)]);
    load_vec(pg3, w2r, &v[VLEN * (icr + ID2)]);
    load_vec(pg3, w2i, &v[VLEN * (ici + ID2)]);

    load_add(pg1, w1r, &v[VLEN * (icr + ID1) + 1]);
    load_add(pg1, w1i, &v[VLEN * (ici + ID1) + 1]);
    load_add(pg1, w2r, &v[VLEN * (icr + ID2) + 1]);
    load_add(pg1, w2i, &v[VLEN * (ici + ID2) + 1]);

    load_vec(pg3, w3r, &v[VLEN * (icr + ID3)]);
    load_vec(pg3, w3i, &v[VLEN * (ici + ID3)]);
    load_vec(pg3, w4r, &v[VLEN * (icr + ID4)]);
    load_vec(pg3, w4i, &v[VLEN * (ici + ID4)]);

    svbool_t pg13 = sveor_z(pg, pg1, pg3);

    load_add(pg1, w3r, &v[VLEN * (icr + ID3) + 1]);
    load_add(pg1, w3i, &v[VLEN * (ici + ID3) + 1]);
    load_add(pg1, w4r, &v[VLEN * (icr + ID4) + 1]);
    load_add(pg1, w4i, &v[VLEN * (ici + ID4) + 1]);

    add_vec(pg13, vt1r, w1r, w4i);
    sub_vec(pg13, vt1i, w1i, w4r);
    add_vec(pg13, vt2r, w2r, w3i);
    sub_vec(pg13, vt2i, w2i, w3r);

    int skip = (VLENY + 1) / 2;
    load_add_gather(pg2, vt1r, &buf[skip * (2 * ic)], idx, skip);
    load_add_gather(pg2, vt1i, &buf[skip * (2 * ic + 1)], idx, skip);
    load_add_gather(pg2, vt2r, &buf[skip * (2 * ic + NVC)], idx, skip);
    load_add_gather(pg2, vt2i, &buf[skip * (2 * ic + 1 + NVC)], idx, skip);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xp2(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3, svuint_t& svidx,
                                 Vsimd_t *v2, REALTYPE *u,
                                 REALTYPE *v1, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp2(pg, pg1, pg2, pg3, vt10, vt11, vt20, vt21, v1, buf, svidx, 0);
    set_sp2_xp2(pg, pg1, pg2, pg3, vt12, vt13, vt22, vt23, v1, buf, svidx, 1);
    set_sp2_xp2(pg, pg1, pg2, pg3, vt14, vt15, vt24, vt25, v1, buf, svidx, 2);

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
  inline void mult_wilson_eo_xp2(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3, svint_t& svidx,
                                 Vsimd_t *v2, REALTYPE *u,
                                 REALTYPE *v1, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp2(pg, pg1, pg2, pg3, vt10, vt11, vt20, vt21, v1, buf, svidx, 0);
    set_sp2_xp2(pg, pg1, pg2, pg3, vt12, vt13, vt22, vt23, v1, buf, svidx, 1);
    set_sp2_xp2(pg, pg1, pg2, pg3, vt14, vt15, vt24, vt25, v1, buf, svidx, 2);

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
  inline void mult_wilson_eo_xp2(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3, svint_t& svidx,
                                 REALTYPE *__restrict v2,
                                 REALTYPE *__restrict u,
                                 REALTYPE *__restrict v1,
                                 REALTYPE *__restrict buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp2(pg, pg1, pg2, pg3, vt10, vt11, vt20, vt21, v1, buf, svidx, 0);
    set_sp2_xp2(pg, pg1, pg2, pg3, vt12, vt13, vt22, vt23, v1, buf, svidx, 1);
    set_sp2_xp2(pg, pg1, pg2, pg3, vt14, vt15, vt24, vt25, v1, buf, svidx, 2);

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
  inline void set_sp2_xp(svbool_t& pg, svbool_t& pg1,
                         svbool_t& pg2, svbool_t& pg3,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *vx, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

    shift_vec_xbw(pg1, pg2, pg3, w1r, &vx[VLEN * (icr + ID1)],
                  &vn[VLEN * (icr + ID1)]);
    shift_vec_xbw(pg1, pg2, pg3, w1i, &vx[VLEN * (ici + ID1)],
                  &vn[VLEN * (ici + ID1)]);

    shift_vec_xbw(pg1, pg2, pg3, w2r, &vx[VLEN * (icr + ID2)],
                  &vn[VLEN * (icr + ID2)]);
    shift_vec_xbw(pg1, pg2, pg3, w2i, &vx[VLEN * (ici + ID2)],
                  &vn[VLEN * (ici + ID2)]);

    shift_vec_xbw(pg1, pg2, pg3, w3r, &vx[VLEN * (icr + ID3)],
                  &vn[VLEN * (icr + ID3)]);
    shift_vec_xbw(pg1, pg2, pg3, w3i, &vx[VLEN * (ici + ID3)],
                  &vn[VLEN * (ici + ID3)]);

    shift_vec_xbw(pg1, pg2, pg3, w4r, &vx[VLEN * (icr + ID4)],
                  &vn[VLEN * (icr + ID4)]);
    shift_vec_xbw(pg1, pg2, pg3, w4i, &vx[VLEN * (ici + ID4)],
                  &vn[VLEN * (ici + ID4)]);
    add_vec(pg, vt1r, w1r, w4i);
    sub_vec(pg, vt1i, w1i, w4r);
    add_vec(pg, vt2r, w2r, w3i);
    sub_vec(pg, vt2i, w2i, w3r);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xpb(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3, Vsimd_t *v2,
                                 REALTYPE *u, REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp(pg, pg1, pg2, pg3, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xp(pg, pg1, pg2, pg3, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xp(pg, pg1, pg2, pg3, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_u(pg, ut10, ut11, ut12, ut13, ut14, ut15, &u[VLEN * (2 * ic)]);
      mult_uv(pg, wt1r, wt1i, ut10, ut11, ut12, ut13, ut14, ut15,
              vt10, vt11, vt12, vt13, vt14, vt15);
      mult_uv(pg, wt2r, wt2i, ut10, ut11, ut12, ut13, ut14, ut15,
              vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_xp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_xp(svbool_t& pg, svbool_t& pg1,
                         svuint_t& idx,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *vx, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

    shift_vec(pg1, idx, w1r, &vx[VLEN * (icr + ID1)],
              &vn[VLEN * (icr + ID1)]);
    shift_vec(pg1, idx, w1i, &vx[VLEN * (ici + ID1)],
              &vn[VLEN * (ici + ID1)]);

    shift_vec(pg1, idx, w2r, &vx[VLEN * (icr + ID2)],
              &vn[VLEN * (icr + ID2)]);
    shift_vec(pg1, idx, w2i, &vx[VLEN * (ici + ID2)],
              &vn[VLEN * (ici + ID2)]);

    shift_vec(pg1, idx, w3r, &vx[VLEN * (icr + ID3)],
              &vn[VLEN * (icr + ID3)]);
    shift_vec(pg1, idx, w3i, &vx[VLEN * (ici + ID3)],
              &vn[VLEN * (ici + ID3)]);

    shift_vec(pg1, idx, w4r, &vx[VLEN * (icr + ID4)],
              &vn[VLEN * (icr + ID4)]);
    shift_vec(pg1, idx, w4i, &vx[VLEN * (ici + ID4)],
              &vn[VLEN * (ici + ID4)]);

    add_vec(pg, vt1r, w1r, w4i);
    sub_vec(pg, vt1i, w1i, w4r);
    add_vec(pg, vt2r, w2r, w3i);
    sub_vec(pg, vt2i, w2i, w3r);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xpb(svbool_t& pg1, svuint_t& idx,
                                 Vsimd_t *v2, REALTYPE *u,
                                 REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp(pg, pg1, idx, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xp(pg, pg1, idx, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xp(pg, pg1, idx, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_u(pg, ut10, ut11, ut12, ut13, ut14, ut15, &u[VLEN * (2 * ic)]);
      mult_uv(pg, wt1r, wt1i, ut10, ut11, ut12, ut13, ut14, ut15,
              vt10, vt11, vt12, vt13, vt14, vt15);
      mult_uv(pg, wt2r, wt2i, ut10, ut11, ut12, ut13, ut14, ut15,
              vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_xp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xpb(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3,
                                 REALTYPE *__restrict v2,
                                 REALTYPE *__restrict u,
                                 REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp(pg, pg1, pg2, pg3, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xp(pg, pg1, pg2, pg3, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xp(pg, pg1, pg2, pg3, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_u(pg, ut10, ut11, ut12, ut13, ut14, ut15, &u[VLEN * (2 * ic)]);
      mult_uv(pg, wt1r, wt1i, ut10, ut11, ut12, ut13, ut14, ut15,
              vt10, vt11, vt12, vt13, vt14, vt15);
      mult_uv(pg, wt2r, wt2i, ut10, ut11, ut12, ut13, ut14, ut15,
              vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_xp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xm1(svbool_t& pg2,
                                 REALTYPE *buf, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm1(pg2, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_xm1(pg2, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_xm1(pg2, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg2, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic]);

      mult_udv(pg2, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg2, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      wt1r = compact_vec(pg2, wt1r);
      wt1i = compact_vec(pg2, wt1i);
      wt2r = compact_vec(pg2, wt2r);
      wt2i = compact_vec(pg2, wt2i);

      int      skip = (VLENY + 1) / 2;
      svbool_t pg1  = set_predicate_whilelt(skip);
      save_vec(pg1, &buf[skip * (2 * ic)], wt1r);
      save_vec(pg1, &buf[skip * (2 * ic + 1)], wt1i);
      save_vec(pg1, &buf[skip * (2 * ic + NVC)], wt2r);
      save_vec(pg1, &buf[skip * (2 * ic + 1 + NVC)], wt2i);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xm1(svbool_t& pg2, svint_t& svidx,
                                 REALTYPE *buf, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm1(pg2, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_xm1(pg2, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_xm1(pg2, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg2, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic]);

      mult_udv(pg2, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg2, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      int skip = (VLENY + 1) / 2;
      save_vec_scatter(pg2, &buf[skip * (2 * ic)], wt1r, svidx);
      save_vec_scatter(pg2, &buf[skip * (2 * ic + 1)], wt1i, svidx);
      save_vec_scatter(pg2, &buf[skip * (2 * ic + NVC)], wt2r, svidx);
      save_vec_scatter(pg2, &buf[skip * (2 * ic + 1 + NVC)], wt2i, svidx);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_xm2(svbool_t& pg1, svbool_t& pg3,
                          svreal_t& vt1r, svreal_t& vt1i,
                          svreal_t& vt2r, svreal_t& vt2i,
                          REALTYPE *vx, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

    load_vec(pg3, w1r, &vx[VLEN * (icr + ID1)]);
    load_add(pg1, w1r, &vx[VLEN * (icr + ID1) - 1]);

    load_vec(pg3, w1i, &vx[VLEN * (ici + ID1)]);
    load_add(pg1, w1i, &vx[VLEN * (ici + ID1) - 1]);

    load_vec(pg3, w2r, &vx[VLEN * (icr + ID2)]);
    load_add(pg1, w2r, &vx[VLEN * (icr + ID2) - 1]);

    load_vec(pg3, w2i, &vx[VLEN * (ici + ID2)]);
    load_add(pg1, w2i, &vx[VLEN * (ici + ID2) - 1]);

    load_vec(pg3, w3r, &vx[VLEN * (icr + ID3)]);
    load_add(pg1, w3r, &vx[VLEN * (icr + ID3) - 1]);

    load_vec(pg3, w3i, &vx[VLEN * (ici + ID3)]);
    load_add(pg1, w3i, &vx[VLEN * (ici + ID3) - 1]);

    load_vec(pg3, w4r, &vx[VLEN * (icr + ID4)]);
    load_add(pg1, w4r, &vx[VLEN * (icr + ID4) - 1]);

    load_vec(pg3, w4i, &vx[VLEN * (ici + ID4)]);
    load_add(pg1, w4i, &vx[VLEN * (ici + ID4) - 1]);

    svbool_t pg0  = set_predicate();
    svbool_t pg13 = sveor_z(pg0, pg1, pg3);
    sub_vec(pg13, vt1r, w1r, w4i);
    add_vec(pg13, vt1i, w1i, w4r);
    sub_vec(pg13, vt2r, w2r, w3i);
    add_vec(pg13, vt2i, w2i, w3r);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xm2(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3, svint_t& svidx,
                                 Vsimd_t *v2, REALTYPE *u,
                                 REALTYPE *v1, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm2(pg1, pg3, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_xm2(pg1, pg3, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_xm2(pg1, pg3, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_xm2_eo(pg1, pg3, ut10, ut11, ut12, ut13, ut14, ut15,
                       &u[VLEN * NVC * ic]);
      svbool_t pg13 = sveor_z(pg, pg1, pg3);
      mult_udv(pg13, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg13, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      int skip = (VLENY + 1) / 2;
      load_add_gather(pg2, wt1r, &buf[skip * (2 * ic)], svidx);
      load_add_gather(pg2, wt1i, &buf[skip * (2 * ic + 1)], svidx);
      load_add_gather(pg2, wt2r, &buf[skip * (2 * ic + NVC)], svidx);
      load_add_gather(pg2, wt2i, &buf[skip * (2 * ic + 1 + NVC)], svidx);

      set_sp4_xm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xm2(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3, svuint_t& idx,
                                 Vsimd_t *v2, REALTYPE *u,
                                 REALTYPE *v1, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm2(pg1, pg3, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_xm2(pg1, pg3, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_xm2(pg1, pg3, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_xm2_eo(pg1, pg3, ut10, ut11, ut12, ut13, ut14, ut15,
                       &u[VLEN * NVC * ic]);
      svbool_t pg13 = sveor_z(pg, pg1, pg3);
      mult_udv(pg13, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg13, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      int skip = (VLENY + 1) / 2;
      load_add_gather(pg2, wt1r, &buf[skip * (2 * ic)], idx, skip);
      load_add_gather(pg2, wt1i, &buf[skip * (2 * ic + 1)], idx, skip);
      load_add_gather(pg2, wt2r, &buf[skip * (2 * ic + NVC)], idx, skip);
      load_add_gather(pg2, wt2i, &buf[skip * (2 * ic + 1 + NVC)], idx, skip);

      set_sp4_xm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xm2(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3, svint_t& svidx,
                                 REALTYPE *__restrict v2,
                                 REALTYPE *__restrict u,
                                 REALTYPE *__restrict v1,
                                 REALTYPE *__restrict buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm2(pg1, pg3, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_xm2(pg1, pg3, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_xm2(pg1, pg3, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_xm2_eo(pg1, pg3, ut10, ut11, ut12, ut13, ut14, ut15,
                       &u[VLEN * NVC * ic]);
      svbool_t pg13 = sveor_z(pg, pg1, pg3);
      mult_udv(pg13, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg13, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      int skip = (VLENY + 1) / 2;
      load_add_gather(pg2, wt1r, &buf[skip * (2 * ic)], svidx);
      load_add_gather(pg2, wt1i, &buf[skip * (2 * ic + 1)], svidx);
      load_add_gather(pg2, wt2r, &buf[skip * (2 * ic + NVC)], svidx);
      load_add_gather(pg2, wt2i, &buf[skip * (2 * ic + 1 + NVC)], svidx);

      set_sp4_xm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_xm(svbool_t& pg, svbool_t& pg1, svbool_t& pg2,
                         svbool_t& pg3,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *vx, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

    shift_vec_xfw(pg1, pg2, pg3, w1r, &vx[VLEN * (icr + ID1)],
                  &vn[VLEN * (icr + ID1)]);
    shift_vec_xfw(pg1, pg2, pg3, w1i, &vx[VLEN * (ici + ID1)],
                  &vn[VLEN * (ici + ID1)]);

    shift_vec_xfw(pg1, pg2, pg3, w2r, &vx[VLEN * (icr + ID2)],
                  &vn[VLEN * (icr + ID2)]);
    shift_vec_xfw(pg1, pg2, pg3, w2i, &vx[VLEN * (ici + ID2)],
                  &vn[VLEN * (ici + ID2)]);

    shift_vec_xfw(pg1, pg2, pg3, w3r, &vx[VLEN * (icr + ID3)],
                  &vn[VLEN * (icr + ID3)]);
    shift_vec_xfw(pg1, pg2, pg3, w3i, &vx[VLEN * (ici + ID3)],
                  &vn[VLEN * (ici + ID3)]);

    shift_vec_xfw(pg1, pg2, pg3, w4r, &vx[VLEN * (icr + ID4)],
                  &vn[VLEN * (icr + ID4)]);
    shift_vec_xfw(pg1, pg2, pg3, w4i, &vx[VLEN * (ici + ID4)],
                  &vn[VLEN * (ici + ID4)]);

    sub_vec(pg, vt1r, w1r, w4i);
    add_vec(pg, vt1i, w1i, w4r);
    sub_vec(pg, vt2r, w2r, w3i);
    add_vec(pg, vt2i, w2i, w3r);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xmb(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3, Vsimd_t *v2,
                                 REALTYPE *u, REALTYPE *un,
                                 REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm(pg, pg1, pg2, pg3, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xm(pg, pg1, pg2, pg3, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xm(pg, pg1, pg2, pg3, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_xm_eo(pg1, pg2, pg3, ut10, ut11, ut12, ut13, ut14, ut15,
                      &u[VLEN * NVC * ic], &un[VLEN * NVC * ic]);

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
  inline void set_sp2_xm(svbool_t& pg, svbool_t& pg1, svuint_t& idx1,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *vx, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

    shift_vec(pg1, idx1, w1r, &vx[VLEN * (icr + ID1)],
              &vn[VLEN * (icr + ID1)]);
    shift_vec(pg1, idx1, w1i, &vx[VLEN * (ici + ID1)],
              &vn[VLEN * (ici + ID1)]);

    shift_vec(pg1, idx1, w2r, &vx[VLEN * (icr + ID2)],
              &vn[VLEN * (icr + ID2)]);
    shift_vec(pg1, idx1, w2i, &vx[VLEN * (ici + ID2)],
              &vn[VLEN * (ici + ID2)]);

    shift_vec(pg1, idx1, w3r, &vx[VLEN * (icr + ID3)],
              &vn[VLEN * (icr + ID3)]);
    shift_vec(pg1, idx1, w3i, &vx[VLEN * (ici + ID3)],
              &vn[VLEN * (ici + ID3)]);

    shift_vec(pg1, idx1, w4r, &vx[VLEN * (icr + ID4)],
              &vn[VLEN * (icr + ID4)]);
    shift_vec(pg1, idx1, w4i, &vx[VLEN * (ici + ID4)],
              &vn[VLEN * (ici + ID4)]);

    sub_vec(pg, vt1r, w1r, w4i);
    add_vec(pg, vt1i, w1i, w4r);
    sub_vec(pg, vt2r, w2r, w3i);
    add_vec(pg, vt2i, w2i, w3r);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_eo_xmb(svbool_t& pg1, svuint_t& idx1,
                                 Vsimd_t *v2,
                                 REALTYPE *u, REALTYPE *un,
                                 REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm(pg, pg1, idx1, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xm(pg, pg1, idx1, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xm(pg, pg1, idx1, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg1, idx1, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic], &un[VLEN * NVC * ic]);

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
  inline void mult_wilson_eo_xmb(svbool_t& pg1, svbool_t& pg2,
                                 svbool_t& pg3,
                                 REALTYPE *__restrict v2,
                                 REALTYPE *u, REALTYPE *un,
                                 REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm(pg, pg1, pg2, pg3, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xm(pg, pg1, pg2, pg3, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xm(pg, pg1, pg2, pg3, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_xm_eo(pg1, pg2, pg3, ut10, ut11, ut12, ut13, ut14, ut15,
                      &u[VLEN * NVC * ic], &un[VLEN * NVC * ic]);

      mult_udv(pg, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);
      set_sp4_xm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }
} // nameless namespace end

#endif
//============================================================END=====
