/*!
      @file    mult_Wilson_parts_qxs_fp16-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2564 $
*/

#ifndef MULT_WILSON_PARTS_QXS_H
#define MULT_WILSON_PARTS_QXS_H

namespace {
//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xp1(svbool_t& pg2, svint_t& svidx,
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

      Vsimd_t y1r, y1i, y2r, y2i;
      save_vec(pg2, &y1r.v[0], v1r);
      save_vec(pg2, &y1i.v[0], v1i);
      save_vec(pg2, &y2r.v[0], v2r);
      save_vec(pg2, &y2i.v[0], v2i);
      save_vec1_x(&buf[VLENY * (2*ic)  ],       &y1r, 0, 1);
      save_vec1_x(&buf[VLENY * (2*ic+1)],       &y1i, 0, 1);
      save_vec1_x(&buf[VLENY * (2*ic   + NVC)], &y2r, 0, 1);
      save_vec1_x(&buf[VLENY * (2*ic+1 + NVC)], &y2i, 0, 1);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_xp2(svbool_t& pg, svbool_t& pg1, svbool_t& pg2,
                          svreal_t& vt1r, svreal_t& vt1i,
                          svreal_t& vt2r, svreal_t& vt2i,
                          REALTYPE *v, REALTYPE *buf, svint_t& index, int ic)
  {
    int icr = ND * 2 * ic;
    int ici = ND * 2 * ic + 1;

    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;
    load_vec(pg1, w1r, &v[VLEN * (icr + ID1) + 1]);
    load_vec(pg1, w1i, &v[VLEN * (ici + ID1) + 1]);
    load_vec(pg1, w2r, &v[VLEN * (icr + ID2) + 1]);
    load_vec(pg1, w2i, &v[VLEN * (ici + ID2) + 1]);
    load_vec(pg1, w3r, &v[VLEN * (icr + ID3) + 1]);
    load_vec(pg1, w3i, &v[VLEN * (ici + ID3) + 1]);
    load_vec(pg1, w4r, &v[VLEN * (icr + ID4) + 1]);
    load_vec(pg1, w4i, &v[VLEN * (ici + ID4) + 1]);

    add_vec(pg1, vt1r, w1r, w4i);
    sub_vec(pg1, vt1i, w1i, w4r);
    add_vec(pg1, vt2r, w2r, w3i);
    sub_vec(pg1, vt2i, w2i, w3r);

    Vsimd_t v1r, v1i, v2r, v2i;
    load_vec1_x(&v1r, VLENX-1, &buf[VLENY * (2*ic)],         1);
    load_vec1_x(&v1i, VLENX-1, &buf[VLENY * (2*ic+1)],       1);
    load_vec1_x(&v2r, VLENX-1, &buf[VLENY * (2*ic   + NVC)], 1);
    load_vec1_x(&v2i, VLENX-1, &buf[VLENY * (2*ic+1 + NVC)], 1);
    load_add(pg2, vt1r, &v1r.v[0]);
    load_add(pg2, vt1i, &v1i.v[0]);
    load_add(pg2, vt2r, &v2r.v[0]);
    load_add(pg2, vt2i, &v2i.v[0]);
  }

//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xp2(svbool_t& pg1, svbool_t& pg2, svint_t& svidx,
                              Vsimd_t *v2, REALTYPE *u,
                              REALTYPE *v1, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp2(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, buf, svidx, 0);
    set_sp2_xp2(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, buf, svidx, 1);
    set_sp2_xp2(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, buf, svidx, 2);

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
  inline void mult_wilson_xp2(svbool_t& pg1, svbool_t& pg2, svint_t& svidx,
                              REALTYPE *__restrict v2,
                              REALTYPE *__restrict u,
                              REALTYPE *__restrict v1,
                              REALTYPE *__restrict buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp2(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, buf, svidx, 0);
    set_sp2_xp2(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, buf, svidx, 1);
    set_sp2_xp2(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, buf, svidx, 2);

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
  inline void set_sp2_xp(svbool_t& pg, svbool_t& pg1, svbool_t& pg2,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *v, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

    shift_vec_xbw(pg1, pg2, w1r, &v[VLEN * (icr + ID1)],
                  &vn[VLEN * (icr + ID1)]);
    shift_vec_xbw(pg1, pg2, w1i, &v[VLEN * (ici + ID1)],
                  &vn[VLEN * (ici + ID1)]);

    shift_vec_xbw(pg1, pg2, w2r, &v[VLEN * (icr + ID2)],
                  &vn[VLEN * (icr + ID2)]);
    shift_vec_xbw(pg1, pg2, w2i, &v[VLEN * (ici + ID2)],
                  &vn[VLEN * (ici + ID2)]);

    shift_vec_xbw(pg1, pg2, w3r, &v[VLEN * (icr + ID3)],
                  &vn[VLEN * (icr + ID3)]);
    shift_vec_xbw(pg1, pg2, w3i, &v[VLEN * (ici + ID3)],
                  &vn[VLEN * (ici + ID3)]);

    shift_vec_xbw(pg1, pg2, w4r, &v[VLEN * (icr + ID4)],
                  &vn[VLEN * (icr + ID4)]);
    shift_vec_xbw(pg1, pg2, w4i, &v[VLEN * (ici + ID4)],
                  &vn[VLEN * (ici + ID4)]);

    add_vec(pg, vt1r, w1r, w4i);
    sub_vec(pg, vt1i, w1i, w4r);
    add_vec(pg, vt2r, w2r, w3i);
    sub_vec(pg, vt2i, w2i, w3r);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xpb(svbool_t& pg1, svbool_t& pg2, Vsimd_t *v2,
                              REALTYPE *u, REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xp(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xp(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);

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
  inline void mult_wilson_xpb(svbool_t& pg1, svbool_t& pg2,
                              REALTYPE *__restrict v2,
                              REALTYPE *u, REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xp(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xp(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xp(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);


    for (int ic = 0; ic < NC; ++ic) {
      svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
      svreal_t wt1r, wt1i, wt2r, wt2i;
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
  inline void set_sp2_xm1(svbool_t& pg,
                          svreal_t& vt1r, svreal_t& vt1i,
                          svreal_t& vt2r, svreal_t& vt2i,
                          REALTYPE *vx, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

    load_vec(pg, w1r, &vx[VLEN * (icr + ID1)]);
    load_vec(pg, w1i, &vx[VLEN * (ici + ID1)]);

    load_vec(pg, w2r, &vx[VLEN * (icr + ID2)]);
    load_vec(pg, w2i, &vx[VLEN * (ici + ID2)]);

    load_vec(pg, w3r, &vx[VLEN * (icr + ID3)]);
    load_vec(pg, w3i, &vx[VLEN * (ici + ID3)]);

    load_vec(pg, w4r, &vx[VLEN * (icr + ID4)]);
    load_vec(pg, w4i, &vx[VLEN * (ici + ID4)]);

    sub_vec(pg, vt1r, w1r, w4i);
    add_vec(pg, vt1i, w1i, w4r);
    sub_vec(pg, vt2r, w2r, w3i);
    add_vec(pg, vt2i, w2i, w3r);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xm1(svbool_t& pg2, svint_t& svidx,
                              REALTYPE *__restrict buf,
                              REALTYPE *__restrict u,
                              REALTYPE *__restrict v1)
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

      Vsimd_t y1r, y1i, y2r, y2i;
      save_vec(pg2, &y1r.v[0], wt1r);
      save_vec(pg2, &y1i.v[0], wt1i);
      save_vec(pg2, &y2r.v[0], wt2r);
      save_vec(pg2, &y2i.v[0], wt2i);
      save_vec1_x(&buf[VLENY * (2*ic)  ],       &y1r, VLENX-1, 1);
      save_vec1_x(&buf[VLENY * (2*ic+1)],       &y1i, VLENX-1, 1);
      save_vec1_x(&buf[VLENY * (2*ic   + NVC)], &y2r, VLENX-1, 1);
      save_vec1_x(&buf[VLENY * (2*ic+1 + NVC)], &y2i, VLENX-1, 1);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_xm2(svbool_t& pg,
                          svreal_t& vt1r, svreal_t& vt1i,
                          svreal_t& vt2r, svreal_t& vt2i,
                          REALTYPE *vx, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

    load_vec(pg, w1r, &vx[VLEN * (icr + ID1) - 1]);
    load_vec(pg, w1i, &vx[VLEN * (ici + ID1) - 1]);
    load_vec(pg, w2r, &vx[VLEN * (icr + ID2) - 1]);
    load_vec(pg, w2i, &vx[VLEN * (ici + ID2) - 1]);
    load_vec(pg, w3r, &vx[VLEN * (icr + ID3) - 1]);
    load_vec(pg, w3i, &vx[VLEN * (ici + ID3) - 1]);
    load_vec(pg, w4r, &vx[VLEN * (icr + ID4) - 1]);
    load_vec(pg, w4i, &vx[VLEN * (ici + ID4) - 1]);

    sub_vec(pg, vt1r, w1r, w4i);
    add_vec(pg, vt1i, w1i, w4r);
    sub_vec(pg, vt2r, w2r, w3i);
    add_vec(pg, vt2i, w2i, w3r);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xm2(svbool_t& pg1, svbool_t& pg2, svint_t& svidx,
                              Vsimd_t *v2, REALTYPE *u,
                              REALTYPE *v1, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm2(pg1, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_xm2(pg1, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_xm2(pg1, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg1, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic - 1]);
      mult_udv(pg1, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg1, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      Vsimd_t v1r, v1i, v2r, v2i;
      load_vec1_x(&v1r, 0, &buf[VLENY * (2*ic)],         1);
      load_vec1_x(&v1i, 0, &buf[VLENY * (2*ic+1)],       1);
      load_vec1_x(&v2r, 0, &buf[VLENY * (2*ic   + NVC)], 1);
      load_vec1_x(&v2i, 0, &buf[VLENY * (2*ic+1 + NVC)], 1);
      load_add(pg2, wt1r, &v1r.v[0]);
      load_add(pg2, wt1i, &v1i.v[0]);
      load_add(pg2, wt2r, &v2r.v[0]);
      load_add(pg2, wt2i, &v2i.v[0]);

      set_sp4_xm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }

//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xm2(svbool_t& pg1, svbool_t& pg2, svint_t& svidx,
                              REALTYPE *v2, REALTYPE *u,
                              REALTYPE *v1, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm2(pg1, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_xm2(pg1, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_xm2(pg1, vt14, vt15, vt24, vt25, v1, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg1, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic - 1]);
      mult_udv(pg1, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg1, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      Vsimd_t v1r, v1i, v2r, v2i;
      load_vec1_x(&v1r, 0, &buf[VLENY * (2*ic)],         1);
      load_vec1_x(&v1i, 0, &buf[VLENY * (2*ic+1)],       1);
      load_vec1_x(&v2r, 0, &buf[VLENY * (2*ic   + NVC)], 1);
      load_vec1_x(&v2i, 0, &buf[VLENY * (2*ic+1 + NVC)], 1);
      load_add(pg2, wt1r, &v1r.v[0]);
      load_add(pg2, wt1i, &v1i.v[0]);
      load_add(pg2, wt2r, &v2r.v[0]);
      load_add(pg2, wt2i, &v2i.v[0]);

      set_sp4_xm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }

//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_xm(svbool_t& pg, svbool_t& pg1, svbool_t& pg2,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *vx, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

    shift_vec_xfw(pg1, pg2, w1r, &vx[VLEN * (icr + ID1)],
                  &vn[VLEN * (icr + ID1)]);
    shift_vec_xfw(pg1, pg2, w1i, &vx[VLEN * (ici + ID1)],
                  &vn[VLEN * (ici + ID1)]);

    shift_vec_xfw(pg1, pg2, w2r, &vx[VLEN * (icr + ID2)],
                  &vn[VLEN * (icr + ID2)]);
    shift_vec_xfw(pg1, pg2, w2i, &vx[VLEN * (ici + ID2)],
                  &vn[VLEN * (ici + ID2)]);

    shift_vec_xfw(pg1, pg2, w3r, &vx[VLEN * (icr + ID3)],
                  &vn[VLEN * (icr + ID3)]);
    shift_vec_xfw(pg1, pg2, w3i, &vx[VLEN * (ici + ID3)],
                  &vn[VLEN * (ici + ID3)]);

    shift_vec_xfw(pg1, pg2, w4r, &vx[VLEN * (icr + ID4)],
                  &vn[VLEN * (icr + ID4)]);
    shift_vec_xfw(pg1, pg2, w4i, &vx[VLEN * (ici + ID4)],
                  &vn[VLEN * (ici + ID4)]);

    sub_vec(pg, vt1r, w1r, w4i);
    add_vec(pg, vt1i, w1i, w4r);
    sub_vec(pg, vt2r, w2r, w3i);
    add_vec(pg, vt2i, w2i, w3r);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_xmb(svbool_t& pg1, svbool_t& pg2, Vsimd_t *v2,
                              REALTYPE *u, REALTYPE *un,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xm(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xm(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_xm(pg1, pg2, ut10, ut11, ut12, ut13, ut14, ut15,
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
  inline void mult_wilson_xmb(svbool_t& pg1, svbool_t& pg2,
                              REALTYPE *__restrict *v2,
                              REALTYPE *u, REALTYPE *un,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xm(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xm(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_xm(pg1, pg2, ut10, ut11, ut12, ut13, ut14, ut15,
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
  inline void mult_wilson_xmb(svbool_t& pg1, svbool_t& pg2, Vsimd_t *v2,
                              Vsimd_t *u,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_xm(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_xm(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_xm(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[NVC * ic].v[0]);

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
  inline void mult_wilson_yp1(svbool_t& pg2,
                              REALTYPE *__restrict buf,
                              REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t vt1r, vt1i, vt2r, vt2i;
      set_sp2_yp(pg2, vt1r, vt1i, vt2r, vt2i, v1, ic);

      save_vec(pg2, &buf[VLENX * (2 * ic)], vt1r);
      save_vec(pg2, &buf[VLENX * (2 * ic + 1)], vt1i);
      save_vec(pg2, &buf[VLENX * (2 * ic + NVC)], vt2r);
      save_vec(pg2, &buf[VLENX * (2 * ic + 1 + NVC)], vt2i);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_yp2(svbool_t& pg1, svbool_t& pg2,
                              Vsimd_t *v2, REALTYPE *u,
                              REALTYPE *v1, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_yp(pg1, vt10, vt11, vt20, vt21, &v1[VLENX], 0);
    set_sp2_yp(pg1, vt12, vt13, vt22, vt23, &v1[VLENX], 1);
    set_sp2_yp(pg1, vt14, vt15, vt24, vt25, &v1[VLENX], 2);

    int offset = -VLENX * (VLENY - 1);
    int ic     = 0;
    load_add(pg2, vt10, &buf[offset + VLENX * (2 * ic)]);
    load_add(pg2, vt11, &buf[offset + VLENX * (2 * ic + 1)]);
    load_add(pg2, vt20, &buf[offset + VLENX * (2 * ic + NVC)]);
    load_add(pg2, vt21, &buf[offset + VLENX * (2 * ic + 1 + NVC)]);
    ic = 1;
    load_add(pg2, vt12, &buf[offset + VLENX * (2 * ic)]);
    load_add(pg2, vt13, &buf[offset + VLENX * (2 * ic + 1)]);
    load_add(pg2, vt22, &buf[offset + VLENX * (2 * ic + NVC)]);
    load_add(pg2, vt23, &buf[offset + VLENX * (2 * ic + 1 + NVC)]);
    ic = 2;
    load_add(pg2, vt14, &buf[offset + VLENX * (2 * ic)]);
    load_add(pg2, vt15, &buf[offset + VLENX * (2 * ic + 1)]);
    load_add(pg2, vt24, &buf[offset + VLENX * (2 * ic + NVC)]);
    load_add(pg2, vt25, &buf[offset + VLENX * (2 * ic + 1 + NVC)]);

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
  inline void mult_wilson_yp2(svbool_t& pg1, svbool_t& pg2,
                              REALTYPE *__restrict v2,
                              REALTYPE *__restrict u,
                              REALTYPE *__restrict v1,
                              REALTYPE *__restrict buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_yp(pg1, vt10, vt11, vt20, vt21, &v1[VLENX], 0);
    set_sp2_yp(pg1, vt12, vt13, vt22, vt23, &v1[VLENX], 1);
    set_sp2_yp(pg1, vt14, vt15, vt24, vt25, &v1[VLENX], 2);

    int offset = -VLENX * (VLENY - 1);
    int ic     = 0;
    load_add(pg2, vt10, &buf[offset + VLENX * (2 * ic)]);
    load_add(pg2, vt11, &buf[offset + VLENX * (2 * ic + 1)]);
    load_add(pg2, vt20, &buf[offset + VLENX * (2 * ic + NVC)]);
    load_add(pg2, vt21, &buf[offset + VLENX * (2 * ic + 1 + NVC)]);
    ic = 1;
    load_add(pg2, vt12, &buf[offset + VLENX * (2 * ic)]);
    load_add(pg2, vt13, &buf[offset + VLENX * (2 * ic + 1)]);
    load_add(pg2, vt22, &buf[offset + VLENX * (2 * ic + NVC)]);
    load_add(pg2, vt23, &buf[offset + VLENX * (2 * ic + 1 + NVC)]);
    ic = 2;
    load_add(pg2, vt14, &buf[offset + VLENX * (2 * ic)]);
    load_add(pg2, vt15, &buf[offset + VLENX * (2 * ic + 1)]);
    load_add(pg2, vt24, &buf[offset + VLENX * (2 * ic + NVC)]);
    load_add(pg2, vt25, &buf[offset + VLENX * (2 * ic + 1 + NVC)]);

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
  inline void set_sp2_yp(svbool_t& pg, svbool_t& pg1, svbool_t& pg2,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *v, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

#if VLENY > 1
    shift_vec_ybw(pg1, pg2, w1r, &v[VLEN * (icr + ID1)],
                  &vn[VLEN * (icr + ID1)]);
    shift_vec_ybw(pg1, pg2, w1i, &v[VLEN * (ici + ID1)],
                  &vn[VLEN * (ici + ID1)]);

    shift_vec_ybw(pg1, pg2, w2r, &v[VLEN * (icr + ID2)],
                  &vn[VLEN * (icr + ID2)]);
    shift_vec_ybw(pg1, pg2, w2i, &v[VLEN * (ici + ID2)],
                  &vn[VLEN * (ici + ID2)]);

    shift_vec_ybw(pg1, pg2, w3r, &v[VLEN * (icr + ID3)],
                  &vn[VLEN * (icr + ID3)]);
    shift_vec_ybw(pg1, pg2, w3i, &v[VLEN * (ici + ID3)],
                  &vn[VLEN * (ici + ID3)]);

    shift_vec_ybw(pg1, pg2, w4r, &v[VLEN * (icr + ID4)],
                  &vn[VLEN * (icr + ID4)]);
    shift_vec_ybw(pg1, pg2, w4i, &v[VLEN * (ici + ID4)],
                  &vn[VLEN * (ici + ID4)]);
#else
    load_vec(pg, w1r, &vn[VLEN * (icr + ID1)]);
    load_vec(pg, w1i, &vn[VLEN * (ici + ID1)]);

    load_vec(pg, w2r, &vn[VLEN * (icr + ID2)]);
    load_vec(pg, w2i, &vn[VLEN * (ici + ID2)]);

    load_vec(pg, w3r, &vn[VLEN * (icr + ID3)]);
    load_vec(pg, w3i, &vn[VLEN * (ici + ID3)]);

    load_vec(pg, w4r, &vn[VLEN * (icr + ID4)]);
    load_vec(pg, w4i, &vn[VLEN * (ici + ID4)]);
#endif

    sub_vec(pg, vt1r, w1r, w4r);
    sub_vec(pg, vt1i, w1i, w4i);
    add_vec(pg, vt2r, w2r, w3r);
    add_vec(pg, vt2i, w2i, w3i);
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_yp(svbool_t& pg, svbool_t& pg1, svuint_t& idx1,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *v, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

#if VLENY > 1
    shift_vec_ybw(w1r, &v[VLEN * (icr + ID1)],
                  &vn[VLEN * (icr + ID1)]);
    shift_vec_ybw(w1i, &v[VLEN * (ici + ID1)],
                  &vn[VLEN * (ici + ID1)]);

    shift_vec_ybw(w2r, &v[VLEN * (icr + ID2)],
                  &vn[VLEN * (icr + ID2)]);
    shift_vec_ybw(w2i, &v[VLEN * (ici + ID2)],
                  &vn[VLEN * (ici + ID2)]);

    shift_vec_ybw(w3r, &v[VLEN * (icr + ID3)],
                  &vn[VLEN * (icr + ID3)]);
    shift_vec_ybw(w3i, &v[VLEN * (ici + ID3)],
                  &vn[VLEN * (ici + ID3)]);

    shift_vec_ybw(w4r, &v[VLEN * (icr + ID4)],
                  &vn[VLEN * (icr + ID4)]);
    shift_vec_ybw(w4i, &v[VLEN * (ici + ID4)],
                  &vn[VLEN * (ici + ID4)]);
#else
    load_vec(pg, w1r, &vn[VLEN * (icr + ID1)]);
    load_vec(pg, w1i, &vn[VLEN * (ici + ID1)]);

    load_vec(pg, w2r, &vn[VLEN * (icr + ID2)]);
    load_vec(pg, w2i, &vn[VLEN * (ici + ID2)]);

    load_vec(pg, w3r, &vn[VLEN * (icr + ID3)]);
    load_vec(pg, w3i, &vn[VLEN * (ici + ID3)]);

    load_vec(pg, w4r, &vn[VLEN * (icr + ID4)]);
    load_vec(pg, w4i, &vn[VLEN * (ici + ID4)]);
#endif

    sub_vec(pg, vt1r, w1r, w4r);
    sub_vec(pg, vt1i, w1i, w4i);
    add_vec(pg, vt2r, w2r, w3r);
    add_vec(pg, vt2i, w2i, w3i);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ypb(svbool_t& pg1, svuint_t& idx1,
                              Vsimd_t *v2, REALTYPE *u,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_yp(pg, pg1, idx1, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_yp(pg, pg1, idx1, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_yp(pg, pg1, idx1, vt14, vt15, vt24, vt25, v1, v1n, 2);

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
  inline void mult_wilson_ypb(svbool_t& pg1, svbool_t& pg2,
                              Vsimd_t *v2, REALTYPE *u,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_yp(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_yp(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_yp(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);

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
  inline void mult_wilson_ypb(svbool_t& pg1, svbool_t& pg2,
                              REALTYPE *__restrict v2,
                              REALTYPE *__restrict u,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    // always v2 != u, but can be v1==v1n
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_yp(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_yp(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_yp(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);


    for (int ic = 0; ic < NC; ++ic) {
      svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
      svreal_t wt1r, wt1i, wt2r, wt2i;

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
  inline void mult_wilson_ym1(svbool_t& pg2,
                              REALTYPE *__restrict buf,
                              REALTYPE *__restrict u,
                              REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;
    set_sp2_ym(pg2, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_ym(pg2, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_ym(pg2, vt14, vt15, vt24, vt25, v1, 2);

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

      int offset = -VLENX * (VLENY - 1);

      save_vec(pg2, &buf[offset + VLENX * (2 * ic)], wt1r);
      save_vec(pg2, &buf[offset + VLENX * (2 * ic + 1)], wt1i);
      save_vec(pg2, &buf[offset + VLENX * (2 * ic + NVC)], wt2r);
      save_vec(pg2, &buf[offset + VLENX * (2 * ic + 1 + NVC)], wt2i);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ym2(svbool_t& pg1, svbool_t& pg2,
                              Vsimd_t *v2, REALTYPE *u,
                              REALTYPE *v1, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_ym(pg1, vt10, vt11, vt20, vt21, &v1[-VLENX], 0);
    set_sp2_ym(pg1, vt12, vt13, vt22, vt23, &v1[-VLENX], 1);
    set_sp2_ym(pg1, vt14, vt15, vt24, vt25, &v1[-VLENX], 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg1, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic - VLENX]);
      mult_udv(pg1, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg1, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      load_add(pg2, wt1r, &buf[VLENX * (2 * ic)]);
      load_add(pg2, wt1i, &buf[VLENX * (2 * ic + 1)]);
      load_add(pg2, wt2r, &buf[VLENX * (2 * ic + NVC)]);
      load_add(pg2, wt2i, &buf[VLENX * (2 * ic + 1 + NVC)]);

      set_sp4_ym(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ym2(svbool_t& pg1, svbool_t& pg2,
                              REALTYPE *__restrict v2,
                              REALTYPE *__restrict u,
                              REALTYPE *__restrict v1,
                              REALTYPE *__restrict buf)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_ym(pg1, vt10, vt11, vt20, vt21, &v1[-VLENX], 0);
    set_sp2_ym(pg1, vt12, vt13, vt22, vt23, &v1[-VLENX], 1);
    set_sp2_ym(pg1, vt14, vt15, vt24, vt25, &v1[-VLENX], 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg1, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[VLEN * NVC * ic - VLENX]);
      mult_udv(pg1, wt1r, wt1i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt10, vt11, vt12, vt13, vt14, vt15);
      mult_udv(pg1, wt2r, wt2i,
               ut10, ut11, ut12, ut13, ut14, ut15,
               vt20, vt21, vt22, vt23, vt24, vt25);

      load_add(pg2, wt1r, &buf[VLENX * (2 * ic)]);
      load_add(pg2, wt1i, &buf[VLENX * (2 * ic + 1)]);
      load_add(pg2, wt2r, &buf[VLENX * (2 * ic + NVC)]);
      load_add(pg2, wt2i, &buf[VLENX * (2 * ic + 1 + NVC)]);

      set_sp4_ym(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_sp2_ym(svbool_t& pg, svbool_t& pg1, svbool_t& pg2,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *vx, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

#if VLENY > 1
    shift_vec_yfw(pg1, pg2, w1r, &vx[VLEN * (icr + ID1)],
                  &vn[VLEN * (icr + ID1)]);
    shift_vec_yfw(pg1, pg2, w1i, &vx[VLEN * (ici + ID1)],
                  &vn[VLEN * (ici + ID1)]);

    shift_vec_yfw(pg1, pg2, w2r, &vx[VLEN * (icr + ID2)],
                  &vn[VLEN * (icr + ID2)]);
    shift_vec_yfw(pg1, pg2, w2i, &vx[VLEN * (ici + ID2)],
                  &vn[VLEN * (ici + ID2)]);

    shift_vec_yfw(pg1, pg2, w3r, &vx[VLEN * (icr + ID3)],
                  &vn[VLEN * (icr + ID3)]);
    shift_vec_yfw(pg1, pg2, w3i, &vx[VLEN * (ici + ID3)],
                  &vn[VLEN * (ici + ID3)]);

    shift_vec_yfw(pg1, pg2, w4r, &vx[VLEN * (icr + ID4)],
                  &vn[VLEN * (icr + ID4)]);
    shift_vec_yfw(pg1, pg2, w4i, &vx[VLEN * (ici + ID4)],
                  &vn[VLEN * (ici + ID4)]);
#else
    load_vec(pg, w1r, &vn[VLEN * (icr + ID1)]);
    load_vec(pg, w1i, &vn[VLEN * (ici + ID1)]);

    load_vec(pg, w2r, &vn[VLEN * (icr + ID2)]);
    load_vec(pg, w2i, &vn[VLEN * (ici + ID2)]);

    load_vec(pg, w3r, &vn[VLEN * (icr + ID3)]);
    load_vec(pg, w3i, &vn[VLEN * (ici + ID3)]);

    load_vec(pg, w4r, &vn[VLEN * (icr + ID4)]);
    load_vec(pg, w4i, &vn[VLEN * (ici + ID4)]);
#endif
    add_vec(pg, vt1r, w1r, w4r);
    add_vec(pg, vt1i, w1i, w4i);
    sub_vec(pg, vt2r, w2r, w3r);
    sub_vec(pg, vt2i, w2i, w3i);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ymb(svbool_t& pg1, svbool_t& pg2,
                              Vsimd_t *v2, REALTYPE *u, REALTYPE *un,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_ym(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_ym(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_ym(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_ym(pg1, pg2, ut10, ut11, ut12, ut13, ut14, ut15,
                   &u[VLEN * NVC * ic], &un[VLEN * NVC * ic]);
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
  template<typename REALTYPE>
  inline void set_sp2_ym(svbool_t& pg, svbool_t& pg1, svuint_t& idx1,
                         svreal_t& vt1r, svreal_t& vt1i,
                         svreal_t& vt2r, svreal_t& vt2i,
                         REALTYPE *vx, REALTYPE *vn, int ic)
  {
    int      icr = ND * 2 * ic;
    int      ici = ND * 2 * ic + 1;
    svreal_t w1r, w1i, w2r, w2i, w3r, w3i, w4r, w4i;

#if VLENY > 1
    shift_vec_yfw(w1r, &vx[VLEN * (icr + ID1)],
                  &vn[VLEN * (icr + ID1)]);
    shift_vec_yfw(w1i, &vx[VLEN * (ici + ID1)],
                  &vn[VLEN * (ici + ID1)]);

    shift_vec_yfw(w2r, &vx[VLEN * (icr + ID2)],
                  &vn[VLEN * (icr + ID2)]);
    shift_vec_yfw(w2i, &vx[VLEN * (ici + ID2)],
                  &vn[VLEN * (ici + ID2)]);

    shift_vec_yfw(w3r, &vx[VLEN * (icr + ID3)],
                  &vn[VLEN * (icr + ID3)]);
    shift_vec_yfw(w3i, &vx[VLEN * (ici + ID3)],
                  &vn[VLEN * (ici + ID3)]);

    shift_vec_yfw(w4r, &vx[VLEN * (icr + ID4)],
                  &vn[VLEN * (icr + ID4)]);
    shift_vec_yfw(w4i, &vx[VLEN * (ici + ID4)],
                  &vn[VLEN * (ici + ID4)]);
#else
    load_vec(pg, w1r, &vn[VLEN * (icr + ID1)]);
    load_vec(pg, w1i, &vn[VLEN * (ici + ID1)]);

    load_vec(pg, w2r, &vn[VLEN * (icr + ID2)]);
    load_vec(pg, w2i, &vn[VLEN * (ici + ID2)]);

    load_vec(pg, w3r, &vn[VLEN * (icr + ID3)]);
    load_vec(pg, w3i, &vn[VLEN * (ici + ID3)]);

    load_vec(pg, w4r, &vn[VLEN * (icr + ID4)]);
    load_vec(pg, w4i, &vn[VLEN * (ici + ID4)]);
#endif
    add_vec(pg, vt1r, w1r, w4r);
    add_vec(pg, vt1i, w1i, w4i);
    sub_vec(pg, vt2r, w2r, w3r);
    sub_vec(pg, vt2i, w2i, w3i);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_ymb(svbool_t& pg1, svuint_t& idx1,
                              Vsimd_t *v2, REALTYPE *u, REALTYPE *un,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_ym(pg, pg1, idx1, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_ym(pg, pg1, idx1, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_ym(pg, pg1, idx1, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_ym(ut10, ut11, ut12, ut13, ut14, ut15,
                   &u[VLEN * NVC * ic], &un[VLEN * NVC * ic]);
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
  template<typename REALTYPE>
  inline void mult_wilson_ymb(svbool_t& pg1, svbool_t& pg2,
                              REALTYPE *__restrict v2,
                              REALTYPE *u, REALTYPE *un,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_ym(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_ym(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_ym(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag_ym(pg1, pg2, ut10, ut11, ut12, ut13, ut14, ut15,
                   &u[VLEN * NVC * ic], &un[VLEN * NVC * ic]);
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
  template<typename REALTYPE>
  inline void mult_wilson_ymb(svbool_t& pg1, svbool_t& pg2,
                              Vsimd_t *v2, Vsimd_t *u,
                              REALTYPE *v1, REALTYPE *v1n)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_ym(pg, pg1, pg2, vt10, vt11, vt20, vt21, v1, v1n, 0);
    set_sp2_ym(pg, pg1, pg2, vt12, vt13, vt22, vt23, v1, v1n, 1);
    set_sp2_ym(pg, pg1, pg2, vt14, vt15, vt24, vt25, v1, v1n, 2);

    svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
    svreal_t wt1r, wt1i, wt2r, wt2i;

    for (int ic = 0; ic < NC; ++ic) {
      load_udag(pg, ut10, ut11, ut12, ut13, ut14, ut15,
                &u[NVC * ic].v[0]);
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
  template<typename REALTYPE>
  inline void mult_wilson_zp1(REALTYPE *__restrict buf,
                              REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t vt1r, vt1i, vt2r, vt2i;
      set_sp2_zp(pg, vt1r, vt1i, vt2r, vt2i, v1, ic);
      save_vec(pg, &buf[VLEN * (2 * ic)], vt1r);
      save_vec(pg, &buf[VLEN * (2 * ic + 1)], vt1i);
      save_vec(pg, &buf[VLEN * (2 * ic + NVC)], vt2r);
      save_vec(pg, &buf[VLEN * (2 * ic + 1 + NVC)], vt2i);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zp2(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf)
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
      set_sp4_zp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zp2(REALTYPE *__restrict v2,
                              REALTYPE *__restrict u,
                              REALTYPE *__restrict buf)
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
      set_sp4_zp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zpb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_zp(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_zp(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_zp(pg, vt14, vt15, vt24, vt25, v1, 2);

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
      set_sp4_zp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zpb(REALTYPE *__restrict v2,
                              REALTYPE *__restrict u,
                              REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_zp(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_zp(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_zp(pg, vt14, vt15, vt24, vt25, v1, 2);

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
      set_sp4_zp(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zm1(REALTYPE *__restrict buf,
                              REALTYPE *__restrict u,
                              REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_zm(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_zm(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_zm(pg, vt14, vt15, vt24, vt25, v1, 2);

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
  inline void mult_wilson_zm2(Vsimd_t *v2, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t wt1r, wt1i, wt2r, wt2i;
      load_vec(pg, wt1r, &buf[VLEN * (2 * ic)]);
      load_vec(pg, wt1i, &buf[VLEN * (2 * ic + 1)]);
      load_vec(pg, wt2r, &buf[VLEN * (2 * ic + NVC)]);
      load_vec(pg, wt2i, &buf[VLEN * (2 * ic + 1 + NVC)]);
      set_sp4_zm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zm2(REALTYPE *__restrict v2,
                              REALTYPE *__restrict buf)
  {
    svbool_t pg = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t wt1r, wt1i, wt2r, wt2i;
      load_vec(pg, wt1r, &buf[VLEN * (2 * ic)]);
      load_vec(pg, wt1i, &buf[VLEN * (2 * ic + 1)]);
      load_vec(pg, wt2r, &buf[VLEN * (2 * ic + NVC)]);
      load_vec(pg, wt2i, &buf[VLEN * (2 * ic + 1 + NVC)]);
      set_sp4_zm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zmb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_zm(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_zm(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_zm(pg, vt14, vt15, vt24, vt25, v1, 2);

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
      set_sp4_zm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_zmb(REALTYPE *__restrict v2,
                              REALTYPE *__restrict u,
                              REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_zm(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_zm(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_zm(pg, vt14, vt15, vt24, vt25, v1, 2);

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
      set_sp4_zm(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tp1_dirac(REALTYPE *__restrict buf,
                                    REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t vt1r, vt1i, vt2r, vt2i;
      set_sp2_tp_dirac(pg, vt1r, vt1i, vt2r, vt2i, v1, ic);
      save_vec(pg, &buf[VLEN * (2 * ic)], vt1r);
      save_vec(pg, &buf[VLEN * (2 * ic + 1)], vt1i);
      save_vec(pg, &buf[VLEN * (2 * ic + NVC)], vt2r);
      save_vec(pg, &buf[VLEN * (2 * ic + 1 + NVC)], vt2i);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tp2_dirac(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf)
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
      set_sp4_tp_dirac(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tp2_dirac(REALTYPE *__restrict v2,
                                    REALTYPE *__restrict u,
                                    REALTYPE *__restrict buf)
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
      set_sp4_tp_dirac(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tpb_dirac(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_tp_dirac(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_tp_dirac(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_tp_dirac(pg, vt14, vt15, vt24, vt25, v1, 2);

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
      set_sp4_tp_dirac(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tpb_dirac(REALTYPE *__restrict v2,
                                    REALTYPE *__restrict u,
                                    REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_tp_dirac(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_tp_dirac(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_tp_dirac(pg, vt14, vt15, vt24, vt25, v1, 2);

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
      set_sp4_tp_dirac(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tm1_dirac(REALTYPE *__restrict buf,
                                    REALTYPE *__restrict u,
                                    REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_tm_dirac(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_tm_dirac(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_tm_dirac(pg, vt14, vt15, vt24, vt25, v1, 2);

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
  inline void mult_wilson_tm2_dirac(Vsimd_t *v2, REALTYPE *buf)
  {
    svbool_t pg = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t wt1r, wt1i, wt2r, wt2i;
      load_vec(pg, wt1r, &buf[VLEN * (2 * ic)]);
      load_vec(pg, wt1i, &buf[VLEN * (2 * ic + 1)]);
      load_vec(pg, wt2r, &buf[VLEN * (2 * ic + NVC)]);
      load_vec(pg, wt2i, &buf[VLEN * (2 * ic + 1 + NVC)]);
      set_sp4_tm_dirac(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tm2_dirac(REALTYPE *__restrict v2,
                                    REALTYPE *__restrict buf)
  {
    svbool_t pg = set_predicate();

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t wt1r, wt1i, wt2r, wt2i;
      load_vec(pg, wt1r, &buf[VLEN * (2 * ic)]);
      load_vec(pg, wt1i, &buf[VLEN * (2 * ic + 1)]);
      load_vec(pg, wt2r, &buf[VLEN * (2 * ic + NVC)]);
      load_vec(pg, wt2i, &buf[VLEN * (2 * ic + 1 + NVC)]);
      set_sp4_tm_dirac(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tmb_dirac(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_tm_dirac(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_tm_dirac(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_tm_dirac(pg, vt14, vt15, vt24, vt25, v1, 2);

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
      set_sp4_tm_dirac(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_tmb_dirac(REALTYPE *__restrict v2,
                                    REALTYPE *__restrict u,
                                    REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
    svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

    set_sp2_tm_dirac(pg, vt10, vt11, vt20, vt21, v1, 0);
    set_sp2_tm_dirac(pg, vt12, vt13, vt22, vt23, v1, 1);
    set_sp2_tm_dirac(pg, vt14, vt15, vt24, vt25, v1, 2);

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
      set_sp4_tm_dirac(pg, v2, wt1r, wt1i, wt2r, wt2i, ic);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_wilson_aypx_save(REALTYPE *__restrict v2, REALTYPE a,
                                    Vsimd_t *__restrict v2v, REALTYPE *__restrict v1)
  {
    svreal_t v2F, v1F;
    svbool_t pg = set_predicate();

    for (int i = 0; i < NVCD; ++i) {
      load_vec(pg, v1F, &v1[VLEN * i]);
      load_vec(pg, v2F, &v2v[i].v[0]);
      // v2F = svmla_m(pg, v1F, v2F, a);  // v1F = v1F + v2F * a
      // save_vec(pg, &v2[VLEN*i], v2F);
      axpy_vec(pg, v1F, a, v2F); // v1F = v1F + v2F * a
      save_vec(pg, &v2[VLEN * i], v1F);
    }
  }


//====================================================================
} // nameless namespace end
#endif
