/*!
      @file    mult_Staggered_parts_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_STAGGERED_PARTS_QXS_INCLUDED
#define MULT_STAGGERED_PARTS_QXS_INCLUDED

//====================================================================
namespace {
  inline void mult_staggered_up(svbool_t pg,
                                svreal_t& vt0, svreal_t& vt1,
                                svreal_t& vt2, svreal_t& vt3,
                                svreal_t& vt4, svreal_t& vt5,
                                real_t *u, real_t *v1)
  {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    load_vec(pg, wt0, &v1[VLEN * 0]);
    load_vec(pg, wt1, &v1[VLEN * 1]);
    load_vec(pg, wt2, &v1[VLEN * 2]);
    load_vec(pg, wt3, &v1[VLEN * 3]);
    load_vec(pg, wt4, &v1[VLEN * 4]);
    load_vec(pg, wt5, &v1[VLEN * 5]);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;
    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 0]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt0, xtr);
    add_vec(pg, vt1, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 2]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt2, xtr);
    add_vec(pg, vt3, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 4]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt4, xtr);
    add_vec(pg, vt5, xti);
  }


  inline void mult_staggered_dn(svbool_t pg,
                                svreal_t& vt0, svreal_t& vt1,
                                svreal_t& vt2, svreal_t& vt3,
                                svreal_t& vt4, svreal_t& vt5,
                                real_t *u, real_t *v1)
  {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    load_vec(pg, wt0, &v1[VLEN * 0]);
    load_vec(pg, wt1, &v1[VLEN * 1]);
    load_vec(pg, wt2, &v1[VLEN * 2]);
    load_vec(pg, wt3, &v1[VLEN * 3]);
    load_vec(pg, wt4, &v1[VLEN * 4]);
    load_vec(pg, wt5, &v1[VLEN * 5]);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;
    load_udag(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 0]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt0, xtr);
    sub_vec(pg, vt1, xti);
    load_udag(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * NVC]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt2, xtr);
    sub_vec(pg, vt3, xti);
    load_udag(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 2 * NVC]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt4, xtr);
    sub_vec(pg, vt5, xti);
  }


  inline void mult_staggered_dn1(svbool_t pg,
                                 real_t *buf, real_t *u, real_t *v1)
  {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    load_vec(pg, wt0, &v1[VLEN * 0]);
    load_vec(pg, wt1, &v1[VLEN * 1]);
    load_vec(pg, wt2, &v1[VLEN * 2]);
    load_vec(pg, wt3, &v1[VLEN * 3]);
    load_vec(pg, wt4, &v1[VLEN * 4]);
    load_vec(pg, wt5, &v1[VLEN * 5]);

    for (int ic = 0; ic < NC; ++ic) {
      svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
      svreal_t xtr, xti;
      load_udag(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * NVC * ic]);
      mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
               wt0, wt1, wt2, wt3, wt4, wt5);
      // svst1(pg, &buf[VLEN * (2*ic)  ], xtr);
      // svst1(pg, &buf[VLEN * (2*ic+1)], xti);
      save_vec(pg, &buf[VLEN * (2 * ic)], xtr);
      save_vec(pg, &buf[VLEN * (2 * ic + 1)], xti);
    }
  }


  inline void mult_staggered_xp(svbool_t pg, svbool_t pg1,
                                svbool_t pg2,
                                svreal_t& vt0, svreal_t& vt1,
                                svreal_t& vt2, svreal_t& vt3,
                                svreal_t& vt4, svreal_t& vt5,
                                real_t *u, real_t *vx, real_t *vn)
  {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    shift_vec_xbw(pg1, pg2, wt0, &vx[VLEN * 0], &vn[VLEN * 0]);
    shift_vec_xbw(pg1, pg2, wt1, &vx[VLEN * 1], &vn[VLEN * 1]);
    shift_vec_xbw(pg1, pg2, wt2, &vx[VLEN * 2], &vn[VLEN * 2]);
    shift_vec_xbw(pg1, pg2, wt3, &vx[VLEN * 3], &vn[VLEN * 3]);
    shift_vec_xbw(pg1, pg2, wt4, &vx[VLEN * 4], &vn[VLEN * 4]);
    shift_vec_xbw(pg1, pg2, wt5, &vx[VLEN * 5], &vn[VLEN * 5]);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;
    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 0]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt0, xtr);
    add_vec(pg, vt1, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 2]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt2, xtr);
    add_vec(pg, vt3, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 4]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt4, xtr);
    add_vec(pg, vt5, xti);
  }


  inline void mult_staggered_xm(svbool_t pg, svbool_t pg1,
                                svbool_t pg2,
                                svreal_t& vt0, svreal_t& vt1,
                                svreal_t& vt2, svreal_t& vt3,
                                svreal_t& vt4, svreal_t& vt5,
                                real_t *ux, real_t *un,
                                real_t *vx, real_t *vn)
  {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    shift_vec_xfw(pg1, pg2, wt0, &vx[VLEN * 0], &vn[VLEN * 0]);
    shift_vec_xfw(pg1, pg2, wt1, &vx[VLEN * 1], &vn[VLEN * 1]);
    shift_vec_xfw(pg1, pg2, wt2, &vx[VLEN * 2], &vn[VLEN * 2]);
    shift_vec_xfw(pg1, pg2, wt3, &vx[VLEN * 3], &vn[VLEN * 3]);
    shift_vec_xfw(pg1, pg2, wt4, &vx[VLEN * 4], &vn[VLEN * 4]);
    shift_vec_xfw(pg1, pg2, wt5, &vx[VLEN * 5], &vn[VLEN * 5]);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;
    load_udag_xm(pg1, pg2, ut0, ut1, ut2, ut3, ut4, ut5,
                 &ux[VLEN * 0], &un[VLEN * 0]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt0, xtr);
    sub_vec(pg, vt1, xti);
    load_udag_xm(pg1, pg2, ut0, ut1, ut2, ut3, ut4, ut5,
                 &ux[VLEN * NVC], &un[VLEN * NVC]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt2, xtr);
    sub_vec(pg, vt3, xti);
    load_udag_xm(pg1, pg2, ut0, ut1, ut2, ut3, ut4, ut5,
                 &ux[VLEN * 2 * NVC], &un[VLEN * 2 * NVC]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt4, xtr);
    sub_vec(pg, vt5, xti);
  }


  inline void mult_staggered_yp(svbool_t pg, svbool_t pg1,
                                svbool_t pg2,
                                svreal_t& vt0, svreal_t& vt1,
                                svreal_t& vt2, svreal_t& vt3,
                                svreal_t& vt4, svreal_t& vt5,
                                real_t *u, real_t *vx, real_t *vn)
  {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    shift_vec_ybw(pg1, pg2, wt0, &vx[VLEN * 0], &vn[VLEN * 0]);
    shift_vec_ybw(pg1, pg2, wt1, &vx[VLEN * 1], &vn[VLEN * 1]);
    shift_vec_ybw(pg1, pg2, wt2, &vx[VLEN * 2], &vn[VLEN * 2]);
    shift_vec_ybw(pg1, pg2, wt3, &vx[VLEN * 3], &vn[VLEN * 3]);
    shift_vec_ybw(pg1, pg2, wt4, &vx[VLEN * 4], &vn[VLEN * 4]);
    shift_vec_ybw(pg1, pg2, wt5, &vx[VLEN * 5], &vn[VLEN * 5]);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;
    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 0]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt0, xtr);
    add_vec(pg, vt1, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 2]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt2, xtr);
    add_vec(pg, vt3, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 4]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt4, xtr);
    add_vec(pg, vt5, xti);
  }


  inline void mult_staggered_ym(svbool_t pg, svbool_t pg1,
                                svbool_t pg2,
                                svreal_t& vt0, svreal_t& vt1,
                                svreal_t& vt2, svreal_t& vt3,
                                svreal_t& vt4, svreal_t& vt5,
                                real_t *ux, real_t *un,
                                real_t *vx, real_t *vn)
  {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    shift_vec_yfw(pg1, pg2, wt0, &vx[VLEN * 0], &vn[VLEN * 0]);
    shift_vec_yfw(pg1, pg2, wt1, &vx[VLEN * 1], &vn[VLEN * 1]);
    shift_vec_yfw(pg1, pg2, wt2, &vx[VLEN * 2], &vn[VLEN * 2]);
    shift_vec_yfw(pg1, pg2, wt3, &vx[VLEN * 3], &vn[VLEN * 3]);
    shift_vec_yfw(pg1, pg2, wt4, &vx[VLEN * 4], &vn[VLEN * 4]);
    shift_vec_yfw(pg1, pg2, wt5, &vx[VLEN * 5], &vn[VLEN * 5]);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;
    load_udag_ym(pg1, pg2, ut0, ut1, ut2, ut3, ut4, ut5,
                 &ux[VLEN * 0], &un[VLEN * 0]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt0, xtr);
    sub_vec(pg, vt1, xti);
    load_udag_ym(pg1, pg2, ut0, ut1, ut2, ut3, ut4, ut5,
                 &ux[VLEN * NVC], &un[VLEN * NVC]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt2, xtr);
    sub_vec(pg, vt3, xti);
    load_udag_ym(pg1, pg2, ut0, ut1, ut2, ut3, ut4, ut5,
                 &ux[VLEN * 2 * NVC], &un[VLEN * 2 * NVC]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt4, xtr);
    sub_vec(pg, vt5, xti);
  }


  inline void mult_staggered_eo_xp(svbool_t pg, svbool_t pg1,
                                   svbool_t pg2, svbool_t pg3,
                                   svreal_t& vt0, svreal_t& vt1,
                                   svreal_t& vt2, svreal_t& vt3,
                                   svreal_t& vt4, svreal_t& vt5,
                                   real_t *u, real_t *vx, real_t *vn)
  {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    shift_vec_xbw(pg1, pg2, pg3, wt0, &vx[VLEN * 0], &vn[VLEN * 0]);
    shift_vec_xbw(pg1, pg2, pg3, wt1, &vx[VLEN * 1], &vn[VLEN * 1]);
    shift_vec_xbw(pg1, pg2, pg3, wt2, &vx[VLEN * 2], &vn[VLEN * 2]);
    shift_vec_xbw(pg1, pg2, pg3, wt3, &vx[VLEN * 3], &vn[VLEN * 3]);
    shift_vec_xbw(pg1, pg2, pg3, wt4, &vx[VLEN * 4], &vn[VLEN * 4]);
    shift_vec_xbw(pg1, pg2, pg3, wt5, &vx[VLEN * 5], &vn[VLEN * 5]);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;
    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 0]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt0, xtr);
    add_vec(pg, vt1, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 2]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt2, xtr);
    add_vec(pg, vt3, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 4]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt4, xtr);
    add_vec(pg, vt5, xti);
  }


  inline void mult_staggered_eo_xm(svbool_t pg, svbool_t pg1,
                                   svbool_t pg2, svbool_t pg3,
                                   svreal_t& vt0, svreal_t& vt1,
                                   svreal_t& vt2, svreal_t& vt3,
                                   svreal_t& vt4, svreal_t& vt5,
                                   real_t *ux, real_t *un,
                                   real_t *vx, real_t *vn)
  {
    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    shift_vec_xfw(pg1, pg2, pg3, wt0, &vx[VLEN * 0], &vn[VLEN * 0]);
    shift_vec_xfw(pg1, pg2, pg3, wt1, &vx[VLEN * 1], &vn[VLEN * 1]);
    shift_vec_xfw(pg1, pg2, pg3, wt2, &vx[VLEN * 2], &vn[VLEN * 2]);
    shift_vec_xfw(pg1, pg2, pg3, wt3, &vx[VLEN * 3], &vn[VLEN * 3]);
    shift_vec_xfw(pg1, pg2, pg3, wt4, &vx[VLEN * 4], &vn[VLEN * 4]);
    shift_vec_xfw(pg1, pg2, pg3, wt5, &vx[VLEN * 5], &vn[VLEN * 5]);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;
    load_udag_xm_eo(pg1, pg2, pg3, ut0, ut1, ut2, ut3, ut4, ut5,
                    &ux[VLEN * 0], &un[VLEN * 0]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt0, xtr);
    sub_vec(pg, vt1, xti);

    load_udag_xm_eo(pg1, pg2, pg3, ut0, ut1, ut2, ut3, ut4, ut5,
                    &ux[VLEN * NVC], &un[VLEN * NVC]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt2, xtr);
    sub_vec(pg, vt3, xti);

    load_udag_xm_eo(pg1, pg2, pg3, ut0, ut1, ut2, ut3, ut4, ut5,
                    &ux[VLEN * 2 * NVC], &un[VLEN * 2 * NVC]);
    mult_udv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    sub_vec(pg, vt4, xtr);
    sub_vec(pg, vt5, xti);
  }


  inline void mult_staggered_eo_xm1(svbool_t pg2, svint_t svidx,
                                    real_t *buf, real_t *u, real_t *v1)
  {
    int Nskipx = (VLENY + 1) / 2;

    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;
    load_vec(pg2, wt0, &v1[VLEN * 0]);
    load_vec(pg2, wt1, &v1[VLEN * 1]);
    load_vec(pg2, wt2, &v1[VLEN * 2]);
    load_vec(pg2, wt3, &v1[VLEN * 3]);
    load_vec(pg2, wt4, &v1[VLEN * 4]);
    load_vec(pg2, wt5, &v1[VLEN * 5]);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;

    load_udag(pg2, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 0]);
    mult_udv(pg2, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    //  svst1_scatter_index(pg2, &buf[Nskipx * 0], svidx, xtr);
    //  svst1_scatter_index(pg2, &buf[Nskipx * 1], svidx, xti);
    save_vec_scatter(pg2, &buf[Nskipx * 0], xtr, svidx);
    save_vec_scatter(pg2, &buf[Nskipx * 1], xti, svidx);

    load_udag(pg2, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * NVC]);
    mult_udv(pg2, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    //svst1_scatter_index(pg2, &buf[Nskipx * 2], svidx, xtr);
    //svst1_scatter_index(pg2, &buf[Nskipx * 3], svidx, xti);
    save_vec_scatter(pg2, &buf[Nskipx * 2], xtr, svidx);
    save_vec_scatter(pg2, &buf[Nskipx * 3], xti, svidx);

    load_udag(pg2, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 2 * NVC]);
    mult_udv(pg2, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    //svst1_scatter_index(pg2, &buf[Nskipx * 4], svidx, xtr);
    //svst1_scatter_index(pg2, &buf[Nskipx * 5], svidx, xti);
    save_vec_scatter(pg2, &buf[Nskipx * 4], xtr, svidx);
    save_vec_scatter(pg2, &buf[Nskipx * 5], xti, svidx);
  }


  inline void mult_staggered_eo_xp2(svbool_t pg1, svbool_t pg2,
                                    svbool_t pg3, svint_t svidx,
                                    svreal_t& vt0, svreal_t& vt1,
                                    svreal_t& vt2, svreal_t& vt3,
                                    svreal_t& vt4, svreal_t& vt5,
                                    real_t *u, real_t *v1, real_t *buf)
  {
    svbool_t pg     = set_predicate();
    int      Nskipx = (VLENY + 1) / 2;

    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;

    load_vec(pg3, wt0, &v1[VLEN * 0]);
    load_add(pg1, wt0, &v1[VLEN * 0 + 1]);
    load_add_gather(pg2, wt0, &buf[Nskipx * 0], svidx);

    load_vec(pg3, wt1, &v1[VLEN * 1]);
    load_add(pg1, wt1, &v1[VLEN * 1 + 1]);
    load_add_gather(pg2, wt1, &buf[Nskipx * 1], svidx);

    load_vec(pg3, wt2, &v1[VLEN * 2]);
    load_add(pg1, wt2, &v1[VLEN * 2 + 1]);
    load_add_gather(pg2, wt2, &buf[Nskipx * 2], svidx);

    load_vec(pg3, wt3, &v1[VLEN * 3]);
    load_add(pg1, wt3, &v1[VLEN * 3 + 1]);
    load_add_gather(pg2, wt3, &buf[Nskipx * 3], svidx);

    load_vec(pg3, wt4, &v1[VLEN * 4]);
    load_add(pg1, wt4, &v1[VLEN * 4 + 1]);
    load_add_gather(pg2, wt4, &buf[Nskipx * 4], svidx);

    load_vec(pg3, wt5, &v1[VLEN * 5]);
    load_add(pg1, wt5, &v1[VLEN * 5 + 1]);
    load_add_gather(pg2, wt5, &buf[Nskipx * 5], svidx);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 0]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt0, xtr);
    add_vec(pg, vt1, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 2]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt2, xtr);
    add_vec(pg, vt3, xti);

    load_u(pg, ut0, ut1, ut2, ut3, ut4, ut5, &u[VLEN * 4]);
    mult_uv(pg, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
            wt0, wt1, wt2, wt3, wt4, wt5);
    add_vec(pg, vt4, xtr);
    add_vec(pg, vt5, xti);
  }


  inline void mult_staggered_eo_xm2(svbool_t pg1, svbool_t pg2,
                                    svbool_t pg3, svint_t svidx,
                                    svreal_t& vt0, svreal_t& vt1,
                                    svreal_t& vt2, svreal_t& vt3,
                                    svreal_t& vt4, svreal_t& vt5,
                                    real_t *u, real_t *v1, real_t *buf)
  {
    svbool_t pg     = set_predicate();
    int      Nskipx = (VLENY + 1) / 2;

    svreal_t wt0, wt1, wt2, wt3, wt4, wt5;

    load_vec(pg3, wt0, &v1[VLEN * 0]);
    load_add(pg1, wt0, &v1[VLEN * 0 - 1]);

    load_vec(pg3, wt1, &v1[VLEN * 1]);
    load_add(pg1, wt1, &v1[VLEN * 1 - 1]);

    load_vec(pg3, wt2, &v1[VLEN * 2]);
    load_add(pg1, wt2, &v1[VLEN * 2 - 1]);

    load_vec(pg3, wt3, &v1[VLEN * 3]);
    load_add(pg1, wt3, &v1[VLEN * 3 - 1]);

    load_vec(pg3, wt4, &v1[VLEN * 4]);
    load_add(pg1, wt4, &v1[VLEN * 4 - 1]);

    load_vec(pg3, wt5, &v1[VLEN * 5]);
    load_add(pg1, wt5, &v1[VLEN * 5 - 1]);

    svbool_t pg13 = sveor_z(pg, pg1, pg3);

    svreal_t ut0, ut1, ut2, ut3, ut4, ut5;
    svreal_t xtr, xti;

    load_udag_xm2_eo(pg1, pg3, ut0, ut1, ut2, ut3, ut4, ut5,
                     &u[VLEN * 0]);
    mult_udv(pg13, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    load_add_gather(pg2, xtr, &buf[Nskipx * 0], svidx);
    load_add_gather(pg2, xti, &buf[Nskipx * 1], svidx);
    sub_vec(pg, vt0, xtr);
    sub_vec(pg, vt1, xti);

    load_udag_xm2_eo(pg1, pg3, ut0, ut1, ut2, ut3, ut4, ut5,
                     &u[VLEN * NVC]);
    mult_udv(pg13, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    load_add_gather(pg2, xtr, &buf[Nskipx * 2], svidx);
    load_add_gather(pg2, xti, &buf[Nskipx * 3], svidx);
    sub_vec(pg, vt2, xtr);
    sub_vec(pg, vt3, xti);

    load_udag_xm2_eo(pg1, pg3, ut0, ut1, ut2, ut3, ut4, ut5,
                     &u[VLEN * 2 * NVC]);
    mult_udv(pg13, xtr, xti, ut0, ut1, ut2, ut3, ut4, ut5,
             wt0, wt1, wt2, wt3, wt4, wt5);
    load_add_gather(pg2, xtr, &buf[Nskipx * 4], svidx);
    load_add_gather(pg2, xti, &buf[Nskipx * 5], svidx);
    sub_vec(pg, vt4, xtr);
    sub_vec(pg, vt5, xti);
  }
} // namespace

#endif
//============================================================END=====
