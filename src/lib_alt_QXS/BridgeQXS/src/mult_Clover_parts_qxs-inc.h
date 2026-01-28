/*!
      @file    mult_Clover_parts_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2514 $
*/

#ifndef MULT_CLOVER_PARTS_QXS_H
#define MULT_CLOVER_PARTS_QXS_H
namespace {
//====================================================================

//11 22 33 44 55 66   12 34 56  23 45   13 24 35 46  15 26  14 36 25  16
  constexpr int clv_idx[72] = {
    0,  -1,  6,  7, 16, 17, 28, 29, 24, 25, 34, 35,
    -1, -1,  1, -1, 12, 13, 18, 19, 32, 33, 26, 27,
    -1, -1, -1, -5,  2, -1,  8,  9, 20, 21, 30, 31,
    -1, -1, -1, -1, -1, -1,  3, -1, 14, 15, 22, 23,
    -1, -1, -1, -1, -1, -1, -1, -1,  4, -1, 10, 11,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  5, -1,
  };

  inline constexpr int index_clv_up(int id1, int ic1, int id2, int ic2, int ri)
  {
    // int i = ic1 + 3 * id1;
    // int j = ic2 + 3 * id2;
    // return clv_idx[2 * (6 * i + j) + ri];
    // N.B. C++ 11 accepts only a return statement for a constexpr function
    // having more statements is a feature of C++14 or later.
    return clv_idx[2 * (6 * (ic1 + 3 * id1) + ic2 + 3 * id2) + ri];
  }


  inline constexpr int index_clv_dn(int id1, int ic1, int id2, int ic2, int ri)
  {
    // int i = ic1 + 3 * id1;
    // int j = ic2 + 3 * id2;
    // return clv_idx[2 * (6 * i + j) + ri] + 36;
    // N.B. C++ 11 accepts only a return statement for a constexpr function
    // having more statements is a feature of C++14 or later.
    return clv_idx[2 * (6 * (ic1 + 3 * id1) + ic2 + 3 * id2) + ri] + 36;
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_clover_csw_aypx(REALTYPE *v2, REALTYPE a, Vsimd_t *v2v,
                                   REALTYPE *ct, REALTYPE *v1)
  {
    svbool_t pg = set_predicate();

    for (int jd = 0; jd < ND2; ++jd) {
      for (int id = 0; id < ND; ++id) {
        int id2 = (id + ND2) % ND;

        int    ig   = VLEN * NDF * (id + ND * jd);
        real_t *ctp = &ct[ig];

        for (int ic = 0; ic < NC; ++ic) {
          int ic2 = 2 * ic;

          svreal_t ut10, ut11, ut12, ut13, ut14, ut15;
          svreal_t vt10, vt11, vt12, vt13, vt14, vt15;
          svreal_t vt20, vt21, vt22, vt23, vt24, vt25;

          load_vec(pg, ut10, &ctp[VLEN * (ic2)]);
          load_vec(pg, ut11, &ctp[VLEN * (ic2 + 1)]);
          load_vec(pg, ut12, &ctp[VLEN * (NVC + ic2)]);
          load_vec(pg, ut13, &ctp[VLEN * (NVC + ic2 + 1)]);
          load_vec(pg, ut14, &ctp[VLEN * (2 * NVC + ic2)]);
          load_vec(pg, ut15, &ctp[VLEN * (2 * NVC + ic2 + 1)]);

          load_vec(pg, vt10, &v1[VLEN * (2 * id)]);
          load_vec(pg, vt11, &v1[VLEN * (2 * id + 1)]);
          load_vec(pg, vt12, &v1[VLEN * (2 * ND + 2 * id)]);
          load_vec(pg, vt13, &v1[VLEN * (2 * ND + 2 * id + 1)]);
          load_vec(pg, vt14, &v1[VLEN * (4 * ND + 2 * id)]);
          load_vec(pg, vt15, &v1[VLEN * (4 * ND + 2 * id + 1)]);

          load_vec(pg, vt20, &v1[VLEN * (2 * id2)]);
          load_vec(pg, vt21, &v1[VLEN * (2 * id2 + 1)]);
          load_vec(pg, vt22, &v1[VLEN * (2 * ND + 2 * id2)]);
          load_vec(pg, vt23, &v1[VLEN * (2 * ND + 2 * id2 + 1)]);
          load_vec(pg, vt24, &v1[VLEN * (4 * ND + 2 * id2)]);
          load_vec(pg, vt25, &v1[VLEN * (4 * ND + 2 * id2 + 1)]);

          svreal_t wt1r, wt1i, wt2r, wt2i;

          mult_uv(pg, wt1r, wt1i,
                  ut10, ut11, ut12, ut13, ut14, ut15,
                  vt10, vt11, vt12, vt13, vt14, vt15);

          mult_uv(pg, wt2r, wt2i,
                  ut10, ut11, ut12, ut13, ut14, ut15,
                  vt20, vt21, vt22, vt23, vt24, vt25);

          int icd1 = 2 * (jd + ND * ic);
          int icd2 = 2 * (jd + ND2 + ND * ic);

          svreal_t xt1r, xt1i;
          load_vec(pg, xt1r, &v2v[icd1].v[0]);
          load_vec(pg, xt1i, &v2v[icd1 + 1].v[0]);
          add_vec(pg, xt1r, wt1r);
          add_vec(pg, xt1i, wt1i);
          save_vec(pg, &v2v[icd1].v[0], xt1r);
          save_vec(pg, &v2v[icd1 + 1].v[0], xt1i);

          svreal_t xt2r, xt2i;
          load_vec(pg, xt2r, &v2v[icd2].v[0]);
          load_vec(pg, xt2i, &v2v[icd2 + 1].v[0]);
          add_vec(pg, xt2r, wt2r);
          add_vec(pg, xt2i, wt2i);
          save_vec(pg, &v2v[icd2].v[0], xt2r);
          save_vec(pg, &v2v[icd2 + 1].v[0], xt2i);
        }
      }
    }

    svreal_t v2F, v1F;
    for (int i = 0; i < NVCD; ++i) {
      load_vec(pg, v1F, &v1[VLEN * i]);
      load_vec(pg, v2F, &v2v[i].v[0]);
      axpy_vec(pg, v1F, a, v2F); // v1F = v1F + v2F * a
      save_vec(pg, &v2[VLEN * i], v1F);
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_2sp_up(svbool_t pg, svreal_t& vtr, svreal_t& vti,
                         REALTYPE *v1, int id, int ic)
  {
    svreal_t wtr, wti;
    load_vec(pg, wtr, &v1[VLEN * (2 * id + 0 + 2 * ND * ic)]);
    load_vec(pg, wti, &v1[VLEN * (2 * id + 1 + 2 * ND * ic)]);
    load_vec(pg, vtr, &v1[VLEN * (2 * id + 4 + 2 * ND * ic)]);
    load_vec(pg, vti, &v1[VLEN * (2 * id + 5 + 2 * ND * ic)]);
    sub_vec(pg, vtr, wtr);
    sub_vec(pg, vti, wti);
  }


//====================================================================
  template<typename REALTYPE>
  inline void set_2sp_dn(svbool_t pg, svreal_t& vtr, svreal_t& vti,
                         REALTYPE *v1, int id, int ic)
  {
    svreal_t wtr, wti;
    load_vec(pg, wtr, &v1[VLEN * (2 * id + 0 + 2 * ND * ic)]);
    load_vec(pg, wti, &v1[VLEN * (2 * id + 1 + 2 * ND * ic)]);
    load_vec(pg, vtr, &v1[VLEN * (2 * id + 4 + 2 * ND * ic)]);
    load_vec(pg, vti, &v1[VLEN * (2 * id + 5 + 2 * ND * ic)]);
    add_vec(pg, vtr, wtr);
    add_vec(pg, vti, wti);
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_clover_csw_aypx_chrot(
    REALTYPE *__restrict v2, REALTYPE a,
    REALTYPE *__restrict v2v,
    REALTYPE *__restrict ct,
    REALTYPE *__restrict v1)
  {
    // v2 = a*(v2v + ct*v1) + v1  [+ chrial <--> dirac ]
    //  std::abort();
    svbool_t pg = set_predicate();

    svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i;
    svreal_t vt4r, vt4i, vt5r, vt5i, vt6r, vt6i;
    { // upper diagonal block
      set_2sp_up(pg, vt1r, vt1i, v1, 0, 0);
      set_2sp_up(pg, vt2r, vt2i, v1, 0, 1);
      set_2sp_up(pg, vt3r, vt3i, v1, 0, 2);
      set_2sp_up(pg, vt4r, vt4i, v1, 1, 0);
      set_2sp_up(pg, vt5r, vt5i, v1, 1, 1);
      set_2sp_up(pg, vt6r, vt6i, v1, 1, 2);

      svreal_t u11r, u22r, u33r, u44r, u55r, u66r;
      load_vec(pg, u11r, &ct[VLEN * index_clv_up(0, 0, 0, 0, 0)]);
      load_vec(pg, u22r, &ct[VLEN * index_clv_up(0, 1, 0, 1, 0)]);
      load_vec(pg, u33r, &ct[VLEN * index_clv_up(0, 2, 0, 2, 0)]);
      load_vec(pg, u44r, &ct[VLEN * index_clv_up(1, 0, 1, 0, 0)]);
      load_vec(pg, u55r, &ct[VLEN * index_clv_up(1, 1, 1, 1, 0)]);
      load_vec(pg, u66r, &ct[VLEN * index_clv_up(1, 2, 1, 2, 0)]);

      svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i;
      svreal_t wt4r, wt4i, wt5r, wt5i, wt6r, wt6i;
      mul_vec(pg, wt1r, u11r, vt1r); //  = U11 x v1
      mul_vec(pg, wt1i, u11r, vt1i);
      mul_vec(pg, wt2r, u22r, vt2r); // += U22 * v2
      mul_vec(pg, wt2i, u22r, vt2i);
      mul_vec(pg, wt3r, u33r, vt3r); // += U33 x v3
      mul_vec(pg, wt3i, u33r, vt3i);

      svreal_t u12r, u12i;
      load_vec(pg, u12r, &ct[VLEN * index_clv_up(0, 0, 0, 1, 0)]);
      load_vec(pg, u12i, &ct[VLEN * index_clv_up(0, 0, 0, 1, 1)]);

      mul_vec(pg, wt4r, u44r, vt4r); // += U44 x v4
      mul_vec(pg, wt4i, u44r, vt4i);
      mul_vec(pg, wt5r, u55r, vt5r); // += U55 x v5
      mul_vec(pg, wt5i, u55r, vt5i);
      mul_vec(pg, wt6r, u66r, vt6r); // += U66 x v6
      mul_vec(pg, wt6i, u66r, vt6i);

      svreal_t u34r, u34i, u56r, u56i;
      load_vec(pg, u34r, &ct[VLEN * index_clv_up(0, 2, 1, 0, 0)]);
      load_vec(pg, u34i, &ct[VLEN * index_clv_up(0, 2, 1, 0, 1)]);
      load_vec(pg, u56r, &ct[VLEN * index_clv_up(1, 1, 1, 2, 0)]);
      load_vec(pg, u56i, &ct[VLEN * index_clv_up(1, 1, 1, 2, 1)]);

      axpy_vec(pg, wt1r, u12r, vt2r); // += U12 x v2
      axpy_vec(pg, wt1i, u12r, vt2i);
      ymax_vec(pg, wt1r, u12i, vt2i);
      axpy_vec(pg, wt1i, u12i, vt2r);
      axpy_vec(pg, wt2r, u12r, vt1r); // = U21 x v1 = (U12)^* x v1
      axpy_vec(pg, wt2i, u12r, vt1i);
      axpy_vec(pg, wt2r, u12i, vt1i);
      ymax_vec(pg, wt2i, u12i, vt1r);

      svreal_t u23r, u23i;
      load_vec(pg, u23r, &ct[VLEN * index_clv_up(0, 1, 0, 2, 0)]);
      load_vec(pg, u23i, &ct[VLEN * index_clv_up(0, 1, 0, 2, 1)]);
      axpy_vec(pg, wt3r, u34r, vt4r); // += U34 x v4
      axpy_vec(pg, wt3i, u34r, vt4i);
      ymax_vec(pg, wt3r, u34i, vt4i);
      axpy_vec(pg, wt3i, u34i, vt4r);
      axpy_vec(pg, wt4r, u34r, vt3r); // += U43 x v3 = (U34)^* x v3
      axpy_vec(pg, wt4i, u34r, vt3i);
      axpy_vec(pg, wt4r, u34i, vt3i);
      ymax_vec(pg, wt4i, u34i, vt3r);

      svreal_t u45r, u45i;
      load_vec(pg, u45r, &ct[VLEN * index_clv_up(1, 0, 1, 1, 0)]);
      load_vec(pg, u45i, &ct[VLEN * index_clv_up(1, 0, 1, 1, 1)]);
      axpy_vec(pg, wt5r, u56r, vt6r); // += U56 x v6
      axpy_vec(pg, wt5i, u56r, vt6i);
      ymax_vec(pg, wt5r, u56i, vt6i);
      axpy_vec(pg, wt5i, u56i, vt6r);
      axpy_vec(pg, wt6r, u56r, vt5r); // += U65 x v5 = (U56)^* x v5
      axpy_vec(pg, wt6i, u56r, vt5i);
      axpy_vec(pg, wt6r, u56i, vt5i);
      ymax_vec(pg, wt6i, u56i, vt5r);

      svreal_t u13r, u13i;
      load_vec(pg, u13r, &ct[VLEN * index_clv_up(0, 0, 0, 2, 0)]);
      load_vec(pg, u13i, &ct[VLEN * index_clv_up(0, 0, 0, 2, 1)]);
      axpy_vec(pg, wt2r, u23r, vt3r); // += U23 x v3
      axpy_vec(pg, wt2i, u23r, vt3i);
      ymax_vec(pg, wt2r, u23i, vt3i);
      axpy_vec(pg, wt2i, u23i, vt3r);
      axpy_vec(pg, wt3r, u23r, vt2r); // += U32 x v2 = (U23)^* x v2
      axpy_vec(pg, wt3i, u23r, vt2i);
      axpy_vec(pg, wt3r, u23i, vt2i);
      ymax_vec(pg, wt3i, u23i, vt2r);

      svreal_t u24r, u24i;
      load_vec(pg, u24r, &ct[VLEN * index_clv_up(0, 1, 1, 0, 0)]);
      load_vec(pg, u24i, &ct[VLEN * index_clv_up(0, 1, 1, 0, 1)]);
      axpy_vec(pg, wt4r, u45r, vt5r); // += U45 x v5
      axpy_vec(pg, wt4i, u45r, vt5i);
      ymax_vec(pg, wt4r, u45i, vt5i);
      axpy_vec(pg, wt4i, u45i, vt5r);
      axpy_vec(pg, wt5r, u45r, vt4r); // += U54 x v4 = (U45)^* x v4
      axpy_vec(pg, wt5i, u45r, vt4i);
      axpy_vec(pg, wt5r, u45i, vt4i);
      ymax_vec(pg, wt5i, u45i, vt4r);

      svreal_t u35r, u35i;
      load_vec(pg, u35r, &ct[VLEN * index_clv_up(0, 2, 1, 1, 0)]);
      load_vec(pg, u35i, &ct[VLEN * index_clv_up(0, 2, 1, 1, 1)]);
      axpy_vec(pg, wt1r, u13r, vt3r); // += U13 x v3
      axpy_vec(pg, wt1i, u13r, vt3i);
      ymax_vec(pg, wt1r, u13i, vt3i);
      axpy_vec(pg, wt1i, u13i, vt3r);
      axpy_vec(pg, wt3r, u13r, vt1r); // = U31 x v1 = (U13)^* x v1
      axpy_vec(pg, wt3i, u13r, vt1i);
      axpy_vec(pg, wt3r, u13i, vt1i);
      ymax_vec(pg, wt3i, u13i, vt1r);

      svreal_t u46r, u46i;
      load_vec(pg, u46r, &ct[VLEN * index_clv_up(1, 0, 1, 2, 0)]);
      load_vec(pg, u46i, &ct[VLEN * index_clv_up(1, 0, 1, 2, 1)]);
      axpy_vec(pg, wt2r, u24r, vt4r); // += U24 x v4
      axpy_vec(pg, wt2i, u24r, vt4i);
      ymax_vec(pg, wt2r, u24i, vt4i);
      axpy_vec(pg, wt2i, u24i, vt4r);
      axpy_vec(pg, wt4r, u24r, vt2r); // += U42 x v2 = (U24)^* x v2
      axpy_vec(pg, wt4i, u24r, vt2i);
      axpy_vec(pg, wt4r, u24i, vt2i);
      ymax_vec(pg, wt4i, u24i, vt2r);

      svreal_t u15r, u15i;
      load_vec(pg, u15r, &ct[VLEN * index_clv_up(0, 0, 1, 1, 0)]);
      load_vec(pg, u15i, &ct[VLEN * index_clv_up(0, 0, 1, 1, 1)]);
      axpy_vec(pg, wt3r, u35r, vt5r); // += U35 x v5
      axpy_vec(pg, wt3i, u35r, vt5i);
      ymax_vec(pg, wt3r, u35i, vt5i);
      axpy_vec(pg, wt3i, u35i, vt5r);
      axpy_vec(pg, wt5r, u35r, vt3r); // += U53 x v3 = (U35)^* x v3
      axpy_vec(pg, wt5i, u35r, vt3i);
      axpy_vec(pg, wt5r, u35i, vt3i);
      ymax_vec(pg, wt5i, u35i, vt3r);

      svreal_t u26r, u26i;
      load_vec(pg, u26r, &ct[VLEN * index_clv_up(0, 1, 1, 2, 0)]);
      load_vec(pg, u26i, &ct[VLEN * index_clv_up(0, 1, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u46r, vt4r); // += U64 x v4 = (U46)^* x v4
      axpy_vec(pg, wt6i, u46r, vt4i);
      axpy_vec(pg, wt6r, u46i, vt4i);
      ymax_vec(pg, wt6i, u46i, vt4r);
      axpy_vec(pg, wt4r, u46r, vt6r); // += U46 x v6
      axpy_vec(pg, wt4i, u46r, vt6i);
      ymax_vec(pg, wt4r, u46i, vt6i);
      axpy_vec(pg, wt4i, u46i, vt6r);

      svreal_t u14r, u14i;
      load_vec(pg, u14r, &ct[VLEN * index_clv_up(0, 0, 1, 0, 0)]);
      load_vec(pg, u14i, &ct[VLEN * index_clv_up(0, 0, 1, 0, 1)]);
      axpy_vec(pg, wt1r, u15r, vt5r); // += U15 x v5
      axpy_vec(pg, wt1i, u15r, vt5i);
      ymax_vec(pg, wt1r, u15i, vt5i);
      axpy_vec(pg, wt1i, u15i, vt5r);
      axpy_vec(pg, wt5r, u15r, vt1r); // = U51 x v1 = (U15)^* x v1
      axpy_vec(pg, wt5i, u15r, vt1i);
      axpy_vec(pg, wt5r, u15i, vt1i);
      ymax_vec(pg, wt5i, u15i, vt1r);

      svreal_t u36r, u36i;
      load_vec(pg, u36r, &ct[VLEN * index_clv_up(0, 2, 1, 2, 0)]);
      load_vec(pg, u36i, &ct[VLEN * index_clv_up(0, 2, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u26r, vt2r); // += U62 x v2 = (U26)^* x v2
      axpy_vec(pg, wt6i, u26r, vt2i);
      axpy_vec(pg, wt6r, u26i, vt2i);
      ymax_vec(pg, wt6i, u26i, vt2r);
      axpy_vec(pg, wt2r, u26r, vt6r); // += U26 x v6
      axpy_vec(pg, wt2i, u26r, vt6i);
      ymax_vec(pg, wt2r, u26i, vt6i);
      axpy_vec(pg, wt2i, u26i, vt6r);

      svreal_t u25r, u25i;
      load_vec(pg, u25r, &ct[VLEN * index_clv_up(0, 1, 1, 1, 0)]);
      load_vec(pg, u25i, &ct[VLEN * index_clv_up(0, 1, 1, 1, 1)]);
      axpy_vec(pg, wt4r, u14r, vt1r); // = U41 x v1 = (U14)^* x v1
      axpy_vec(pg, wt4i, u14r, vt1i);
      axpy_vec(pg, wt4r, u14i, vt1i);
      ymax_vec(pg, wt4i, u14i, vt1r);
      axpy_vec(pg, wt1r, u14r, vt4r); // += U14 x v4
      axpy_vec(pg, wt1i, u14r, vt4i);
      ymax_vec(pg, wt1r, u14i, vt4i);
      axpy_vec(pg, wt1i, u14i, vt4r);

      svreal_t u16r, u16i;
      load_vec(pg, u16r, &ct[VLEN * index_clv_up(0, 0, 1, 2, 0)]);
      load_vec(pg, u16i, &ct[VLEN * index_clv_up(0, 0, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u36r, vt3r); // += U63 x v3 = (U36)^* x v3
      axpy_vec(pg, wt6i, u36r, vt3i);
      axpy_vec(pg, wt6r, u36i, vt3i);
      ymax_vec(pg, wt6i, u36i, vt3r);
      axpy_vec(pg, wt3r, u36r, vt6r); // += U36 x v6
      axpy_vec(pg, wt3i, u36r, vt6i);
      ymax_vec(pg, wt3r, u36i, vt6i);
      axpy_vec(pg, wt3i, u36i, vt6r);


      axpy_vec(pg, wt2r, u25r, vt5r); // += U25 x v5
      axpy_vec(pg, wt2i, u25r, vt5i);
      ymax_vec(pg, wt2r, u25i, vt5i);
      axpy_vec(pg, wt2i, u25i, vt5r);
      axpy_vec(pg, wt5r, u25r, vt2r); // += U52 x v2 = (U25)^* x v2
      axpy_vec(pg, wt5i, u25r, vt2i);
      axpy_vec(pg, wt5r, u25i, vt2i);
      ymax_vec(pg, wt5i, u25i, vt2r);

      axpy_vec(pg, wt1r, u16r, vt6r); // += U16 x v6
      axpy_vec(pg, wt1i, u16r, vt6i);
      ymax_vec(pg, wt1r, u16i, vt6i);
      axpy_vec(pg, wt1i, u16i, vt6r);
      axpy_vec(pg, wt6r, u16r, vt1r); // = U61 x v1 = (U16)^* x v1
      axpy_vec(pg, wt6i, u16r, vt1i);
      axpy_vec(pg, wt6r, u16i, vt1i);
      ymax_vec(pg, wt6i, u16i, vt1r);

      {
        svreal_t yt1r, yt1i, yt2r, yt2i, yt3r, yt3i;
        load_vec(pg, yt1r, &v2v[VLEN * (2 * (2 + 4 * 0))]);
        load_vec(pg, yt1i, &v2v[VLEN * (1 + 2 * (2 + 4 * 0))]);
        load_vec(pg, yt2r, &v2v[VLEN * (2 * (2 + 4 * 1))]);
        load_vec(pg, yt2i, &v2v[VLEN * (1 + 2 * (2 + 4 * 1))]);
        load_vec(pg, yt3r, &v2v[VLEN * (2 * (2 + 4 * 2))]);
        load_vec(pg, yt3i, &v2v[VLEN * (1 + 2 * (2 + 4 * 2))]);
        add_vec(pg, yt1r, wt1r);
        add_vec(pg, yt1i, wt1i);
        add_vec(pg, yt2r, wt2r);
        add_vec(pg, yt2i, wt2i);
        add_vec(pg, yt3r, wt3r);
        add_vec(pg, yt3i, wt3i);
        save_vec(pg, &v2v[VLEN * (2 * (2 + 4 * 0))], yt1r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (2 + 4 * 0))], yt1i);
        save_vec(pg, &v2v[VLEN * (2 * (2 + 4 * 1))], yt2r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (2 + 4 * 1))], yt2i);
        save_vec(pg, &v2v[VLEN * (2 * (2 + 4 * 2))], yt3r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (2 + 4 * 2))], yt3i);

        svreal_t yt4r, yt4i, yt5r, yt5i, yt6r, yt6i;
        load_vec(pg, yt4r, &v2v[VLEN * (2 * (3 + 4 * 0))]);
        load_vec(pg, yt4i, &v2v[VLEN * (1 + 2 * (3 + 4 * 0))]);
        load_vec(pg, yt5r, &v2v[VLEN * (2 * (3 + 4 * 1))]);
        load_vec(pg, yt5i, &v2v[VLEN * (1 + 2 * (3 + 4 * 1))]);
        load_vec(pg, yt6r, &v2v[VLEN * (2 * (3 + 4 * 2))]);
        load_vec(pg, yt6i, &v2v[VLEN * (1 + 2 * (3 + 4 * 2))]);
        add_vec(pg, yt4r, wt4r);
        add_vec(pg, yt4i, wt4i);
        add_vec(pg, yt5r, wt5r);
        add_vec(pg, yt5i, wt5i);
        add_vec(pg, yt6r, wt6r);
        add_vec(pg, yt6i, wt6i);
        save_vec(pg, &v2v[VLEN * (2 * (3 + 4 * 0))], yt4r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (3 + 4 * 0))], yt4i);
        save_vec(pg, &v2v[VLEN * (2 * (3 + 4 * 1))], yt5r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (3 + 4 * 1))], yt5i);
        save_vec(pg, &v2v[VLEN * (2 * (3 + 4 * 2))], yt6r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (3 + 4 * 2))], yt6i);
      }
      {
        svreal_t yt1r, yt1i, yt2r, yt2i, yt3r, yt3i;
        load_vec(pg, yt1r, &v2v[VLEN * (2 * (0 + 4 * 0))]);
        load_vec(pg, yt1i, &v2v[VLEN * (1 + 2 * (0 + 4 * 0))]);
        load_vec(pg, yt2r, &v2v[VLEN * (2 * (0 + 4 * 1))]);
        load_vec(pg, yt2i, &v2v[VLEN * (1 + 2 * (0 + 4 * 1))]);
        load_vec(pg, yt3r, &v2v[VLEN * (2 * (0 + 4 * 2))]);
        load_vec(pg, yt3i, &v2v[VLEN * (1 + 2 * (0 + 4 * 2))]);
        sub_vec(pg, yt1r, wt1r);
        sub_vec(pg, yt1i, wt1i);
        sub_vec(pg, yt2r, wt2r);
        sub_vec(pg, yt2i, wt2i);
        sub_vec(pg, yt3r, wt3r);
        sub_vec(pg, yt3i, wt3i);
        save_vec(pg, &v2v[VLEN * (2 * (0 + 4 * 0))], yt1r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (0 + 4 * 0))], yt1i);
        save_vec(pg, &v2v[VLEN * (2 * (0 + 4 * 1))], yt2r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (0 + 4 * 1))], yt2i);
        save_vec(pg, &v2v[VLEN * (2 * (0 + 4 * 2))], yt3r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (0 + 4 * 2))], yt3i);

        svreal_t yt4r, yt4i, yt5r, yt5i, yt6r, yt6i;
        load_vec(pg, yt4r, &v2v[VLEN * (2 * (1 + 4 * 0))]);
        load_vec(pg, yt4i, &v2v[VLEN * (1 + 2 * (1 + 4 * 0))]);
        load_vec(pg, yt5r, &v2v[VLEN * (2 * (1 + 4 * 1))]);
        load_vec(pg, yt5i, &v2v[VLEN * (1 + 2 * (1 + 4 * 1))]);
        load_vec(pg, yt6r, &v2v[VLEN * (2 * (1 + 4 * 2))]);
        load_vec(pg, yt6i, &v2v[VLEN * (1 + 2 * (1 + 4 * 2))]);
        sub_vec(pg, yt4r, wt4r);
        sub_vec(pg, yt4i, wt4i);
        sub_vec(pg, yt5r, wt5r);
        sub_vec(pg, yt5i, wt5i);
        sub_vec(pg, yt6r, wt6r);
        sub_vec(pg, yt6i, wt6i);
        save_vec(pg, &v2v[VLEN * (2 * (1 + 4 * 0))], yt4r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (1 + 4 * 0))], yt4i);
        save_vec(pg, &v2v[VLEN * (2 * (1 + 4 * 1))], yt5r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (1 + 4 * 1))], yt5i);
        save_vec(pg, &v2v[VLEN * (2 * (1 + 4 * 2))], yt6r);
        save_vec(pg, &v2v[VLEN * (1 + 2 * (1 + 4 * 2))], yt6i);
      }
    }
    { // lower diagonal block
      set_2sp_dn(pg, vt1r, vt1i, v1, 0, 0);
      set_2sp_dn(pg, vt2r, vt2i, v1, 0, 1);
      set_2sp_dn(pg, vt3r, vt3i, v1, 0, 2);
      set_2sp_dn(pg, vt4r, vt4i, v1, 1, 0);
      set_2sp_dn(pg, vt5r, vt5i, v1, 1, 1);
      set_2sp_dn(pg, vt6r, vt6i, v1, 1, 2);

      svreal_t u11r, u22r, u33r, u44r, u55r, u66r;
      load_vec(pg, u11r, &ct[VLEN * index_clv_dn(0, 0, 0, 0, 0)]);
      load_vec(pg, u22r, &ct[VLEN * index_clv_dn(0, 1, 0, 1, 0)]);
      load_vec(pg, u33r, &ct[VLEN * index_clv_dn(0, 2, 0, 2, 0)]);
      load_vec(pg, u44r, &ct[VLEN * index_clv_dn(1, 0, 1, 0, 0)]);
      load_vec(pg, u55r, &ct[VLEN * index_clv_dn(1, 1, 1, 1, 0)]);
      load_vec(pg, u66r, &ct[VLEN * index_clv_dn(1, 2, 1, 2, 0)]);



      svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i;
      svreal_t wt4r, wt4i, wt5r, wt5i, wt6r, wt6i;
      mul_vec(pg, wt1r, u11r, vt1r); //  = U11 x v1
      mul_vec(pg, wt1i, u11r, vt1i);
      mul_vec(pg, wt2r, u22r, vt2r); // += U22 * v2
      mul_vec(pg, wt2i, u22r, vt2i);
      mul_vec(pg, wt3r, u33r, vt3r); // += U33 x v3
      mul_vec(pg, wt3i, u33r, vt3i);

      svreal_t u12r, u12i;
      load_vec(pg, u12r, &ct[VLEN * index_clv_dn(0, 0, 0, 1, 0)]);
      load_vec(pg, u12i, &ct[VLEN * index_clv_dn(0, 0, 0, 1, 1)]);

      mul_vec(pg, wt4r, u44r, vt4r); // += U44 x v4
      mul_vec(pg, wt4i, u44r, vt4i);
      mul_vec(pg, wt5r, u55r, vt5r); // += U55 x v5
      mul_vec(pg, wt5i, u55r, vt5i);
      mul_vec(pg, wt6r, u66r, vt6r); // += U66 x v6
      mul_vec(pg, wt6i, u66r, vt6i);

      svreal_t u34r, u34i, u56r, u56i;
      load_vec(pg, u34r, &ct[VLEN * index_clv_dn(0, 2, 1, 0, 0)]);
      load_vec(pg, u34i, &ct[VLEN * index_clv_dn(0, 2, 1, 0, 1)]);
      load_vec(pg, u56r, &ct[VLEN * index_clv_dn(1, 1, 1, 2, 0)]);
      load_vec(pg, u56i, &ct[VLEN * index_clv_dn(1, 1, 1, 2, 1)]);

      axpy_vec(pg, wt1r, u12r, vt2r); // += U12 x v2
      axpy_vec(pg, wt1i, u12r, vt2i);
      ymax_vec(pg, wt1r, u12i, vt2i);
      axpy_vec(pg, wt1i, u12i, vt2r);
      axpy_vec(pg, wt2r, u12r, vt1r); // = U21 x v1 = (U12)^* x v1
      axpy_vec(pg, wt2i, u12r, vt1i);
      axpy_vec(pg, wt2r, u12i, vt1i);
      ymax_vec(pg, wt2i, u12i, vt1r);

      svreal_t u23r, u23i;
      load_vec(pg, u23r, &ct[VLEN * index_clv_dn(0, 1, 0, 2, 0)]);
      load_vec(pg, u23i, &ct[VLEN * index_clv_dn(0, 1, 0, 2, 1)]);
      axpy_vec(pg, wt3r, u34r, vt4r); // += U34 x v4
      axpy_vec(pg, wt3i, u34r, vt4i);
      ymax_vec(pg, wt3r, u34i, vt4i);
      axpy_vec(pg, wt3i, u34i, vt4r);
      axpy_vec(pg, wt4r, u34r, vt3r); // += U43 x v3 = (U34)^* x v3
      axpy_vec(pg, wt4i, u34r, vt3i);
      axpy_vec(pg, wt4r, u34i, vt3i);
      ymax_vec(pg, wt4i, u34i, vt3r);

      svreal_t u45r, u45i;
      load_vec(pg, u45r, &ct[VLEN * index_clv_dn(1, 0, 1, 1, 0)]);
      load_vec(pg, u45i, &ct[VLEN * index_clv_dn(1, 0, 1, 1, 1)]);
      axpy_vec(pg, wt5r, u56r, vt6r); // += U56 x v6
      axpy_vec(pg, wt5i, u56r, vt6i);
      ymax_vec(pg, wt5r, u56i, vt6i);
      axpy_vec(pg, wt5i, u56i, vt6r);
      axpy_vec(pg, wt6r, u56r, vt5r); // += U65 x v5 = (U56)^* x v5
      axpy_vec(pg, wt6i, u56r, vt5i);
      axpy_vec(pg, wt6r, u56i, vt5i);
      ymax_vec(pg, wt6i, u56i, vt5r);

      svreal_t u13r, u13i;
      load_vec(pg, u13r, &ct[VLEN * index_clv_dn(0, 0, 0, 2, 0)]);
      load_vec(pg, u13i, &ct[VLEN * index_clv_dn(0, 0, 0, 2, 1)]);
      axpy_vec(pg, wt2r, u23r, vt3r); // += U23 x v3
      axpy_vec(pg, wt2i, u23r, vt3i);
      ymax_vec(pg, wt2r, u23i, vt3i);
      axpy_vec(pg, wt2i, u23i, vt3r);
      axpy_vec(pg, wt3r, u23r, vt2r); // += U32 x v2 = (U23)^* x v2
      axpy_vec(pg, wt3i, u23r, vt2i);
      axpy_vec(pg, wt3r, u23i, vt2i);
      ymax_vec(pg, wt3i, u23i, vt2r);

      svreal_t u24r, u24i;
      load_vec(pg, u24r, &ct[VLEN * index_clv_dn(0, 1, 1, 0, 0)]);
      load_vec(pg, u24i, &ct[VLEN * index_clv_dn(0, 1, 1, 0, 1)]);
      axpy_vec(pg, wt4r, u45r, vt5r); // += U45 x v5
      axpy_vec(pg, wt4i, u45r, vt5i);
      ymax_vec(pg, wt4r, u45i, vt5i);
      axpy_vec(pg, wt4i, u45i, vt5r);
      axpy_vec(pg, wt5r, u45r, vt4r); // += U54 x v4 = (U45)^* x v4
      axpy_vec(pg, wt5i, u45r, vt4i);
      axpy_vec(pg, wt5r, u45i, vt4i);
      ymax_vec(pg, wt5i, u45i, vt4r);

      svreal_t u35r, u35i;
      load_vec(pg, u35r, &ct[VLEN * index_clv_dn(0, 2, 1, 1, 0)]);
      load_vec(pg, u35i, &ct[VLEN * index_clv_dn(0, 2, 1, 1, 1)]);
      axpy_vec(pg, wt1r, u13r, vt3r); // += U13 x v3
      axpy_vec(pg, wt1i, u13r, vt3i);
      ymax_vec(pg, wt1r, u13i, vt3i);
      axpy_vec(pg, wt1i, u13i, vt3r);
      axpy_vec(pg, wt3r, u13r, vt1r); // = U31 x v1 = (U13)^* x v1
      axpy_vec(pg, wt3i, u13r, vt1i);
      axpy_vec(pg, wt3r, u13i, vt1i);
      ymax_vec(pg, wt3i, u13i, vt1r);

      svreal_t u46r, u46i;
      load_vec(pg, u46r, &ct[VLEN * index_clv_dn(1, 0, 1, 2, 0)]);
      load_vec(pg, u46i, &ct[VLEN * index_clv_dn(1, 0, 1, 2, 1)]);
      axpy_vec(pg, wt2r, u24r, vt4r); // += U24 x v4
      axpy_vec(pg, wt2i, u24r, vt4i);
      ymax_vec(pg, wt2r, u24i, vt4i);
      axpy_vec(pg, wt2i, u24i, vt4r);
      axpy_vec(pg, wt4r, u24r, vt2r); // += U42 x v2 = (U24)^* x v2
      axpy_vec(pg, wt4i, u24r, vt2i);
      axpy_vec(pg, wt4r, u24i, vt2i);
      ymax_vec(pg, wt4i, u24i, vt2r);

      svreal_t u15r, u15i;
      load_vec(pg, u15r, &ct[VLEN * index_clv_dn(0, 0, 1, 1, 0)]);
      load_vec(pg, u15i, &ct[VLEN * index_clv_dn(0, 0, 1, 1, 1)]);
      axpy_vec(pg, wt3r, u35r, vt5r); // += U35 x v5
      axpy_vec(pg, wt3i, u35r, vt5i);
      ymax_vec(pg, wt3r, u35i, vt5i);
      axpy_vec(pg, wt3i, u35i, vt5r);
      axpy_vec(pg, wt5r, u35r, vt3r); // += U53 x v3 = (U35)^* x v3
      axpy_vec(pg, wt5i, u35r, vt3i);
      axpy_vec(pg, wt5r, u35i, vt3i);
      ymax_vec(pg, wt5i, u35i, vt3r);

      svreal_t u26r, u26i;
      load_vec(pg, u26r, &ct[VLEN * index_clv_dn(0, 1, 1, 2, 0)]);
      load_vec(pg, u26i, &ct[VLEN * index_clv_dn(0, 1, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u46r, vt4r); // += U64 x v4 = (U46)^* x v4
      axpy_vec(pg, wt6i, u46r, vt4i);
      axpy_vec(pg, wt6r, u46i, vt4i);
      ymax_vec(pg, wt6i, u46i, vt4r);
      axpy_vec(pg, wt4r, u46r, vt6r); // += U46 x v6
      axpy_vec(pg, wt4i, u46r, vt6i);
      ymax_vec(pg, wt4r, u46i, vt6i);
      axpy_vec(pg, wt4i, u46i, vt6r);

      svreal_t u14r, u14i;
      load_vec(pg, u14r, &ct[VLEN * index_clv_dn(0, 0, 1, 0, 0)]);
      load_vec(pg, u14i, &ct[VLEN * index_clv_dn(0, 0, 1, 0, 1)]);
      axpy_vec(pg, wt1r, u15r, vt5r); // += U15 x v5
      axpy_vec(pg, wt1i, u15r, vt5i);
      ymax_vec(pg, wt1r, u15i, vt5i);
      axpy_vec(pg, wt1i, u15i, vt5r);
      axpy_vec(pg, wt5r, u15r, vt1r); // = U51 x v1 = (U15)^* x v1
      axpy_vec(pg, wt5i, u15r, vt1i);
      axpy_vec(pg, wt5r, u15i, vt1i);
      ymax_vec(pg, wt5i, u15i, vt1r);

      svreal_t u36r, u36i;
      load_vec(pg, u36r, &ct[VLEN * index_clv_dn(0, 2, 1, 2, 0)]);
      load_vec(pg, u36i, &ct[VLEN * index_clv_dn(0, 2, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u26r, vt2r); // += U62 x v2 = (U26)^* x v2
      axpy_vec(pg, wt6i, u26r, vt2i);
      axpy_vec(pg, wt6r, u26i, vt2i);
      ymax_vec(pg, wt6i, u26i, vt2r);
      axpy_vec(pg, wt2r, u26r, vt6r); // += U26 x v6
      axpy_vec(pg, wt2i, u26r, vt6i);
      ymax_vec(pg, wt2r, u26i, vt6i);
      axpy_vec(pg, wt2i, u26i, vt6r);

      svreal_t u25r, u25i;
      load_vec(pg, u25r, &ct[VLEN * index_clv_dn(0, 1, 1, 1, 0)]);
      load_vec(pg, u25i, &ct[VLEN * index_clv_dn(0, 1, 1, 1, 1)]);
      axpy_vec(pg, wt4r, u14r, vt1r); // = U41 x v1 = (U14)^* x v1
      axpy_vec(pg, wt4i, u14r, vt1i);
      axpy_vec(pg, wt4r, u14i, vt1i);
      ymax_vec(pg, wt4i, u14i, vt1r);
      axpy_vec(pg, wt1r, u14r, vt4r); // += U14 x v4
      axpy_vec(pg, wt1i, u14r, vt4i);
      ymax_vec(pg, wt1r, u14i, vt4i);
      axpy_vec(pg, wt1i, u14i, vt4r);

      svreal_t u16r, u16i;
      load_vec(pg, u16r, &ct[VLEN * index_clv_dn(0, 0, 1, 2, 0)]);
      load_vec(pg, u16i, &ct[VLEN * index_clv_dn(0, 0, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u36r, vt3r); // += U63 x v3 = (U36)^* x v3
      axpy_vec(pg, wt6i, u36r, vt3i);
      axpy_vec(pg, wt6r, u36i, vt3i);
      ymax_vec(pg, wt6i, u36i, vt3r);
      axpy_vec(pg, wt3r, u36r, vt6r); // += U36 x v6
      axpy_vec(pg, wt3i, u36r, vt6i);
      ymax_vec(pg, wt3r, u36i, vt6i);
      axpy_vec(pg, wt3i, u36i, vt6r);


      axpy_vec(pg, wt2r, u25r, vt5r); // += U25 x v5
      axpy_vec(pg, wt2i, u25r, vt5i);
      ymax_vec(pg, wt2r, u25i, vt5i);
      axpy_vec(pg, wt2i, u25i, vt5r);
      axpy_vec(pg, wt5r, u25r, vt2r); // += U52 x v2 = (U25)^* x v2
      axpy_vec(pg, wt5i, u25r, vt2i);
      axpy_vec(pg, wt5r, u25i, vt2i);
      ymax_vec(pg, wt5i, u25i, vt2r);

      axpy_vec(pg, wt1r, u16r, vt6r); // += U16 x v6
      axpy_vec(pg, wt1i, u16r, vt6i);
      ymax_vec(pg, wt1r, u16i, vt6i);
      axpy_vec(pg, wt1i, u16i, vt6r);
      axpy_vec(pg, wt6r, u16r, vt1r); // = U61 x v1 = (U16)^* x v1
      axpy_vec(pg, wt6i, u16r, vt1i);
      axpy_vec(pg, wt6r, u16i, vt1i);
      ymax_vec(pg, wt6i, u16i, vt1r);


      for (int sp = 0; sp < 4; sp += 2) {
        { // store to sp1, sp3
          svreal_t yt1r, yt1i, yt2r, yt2i, yt3r, yt3i;
          load_vec(pg, yt1r, &v2v[VLEN * (2 * (sp + 4 * 0))]);
          load_vec(pg, yt1i, &v2v[VLEN * (1 + 2 * (sp + 4 * 0))]);
          load_vec(pg, yt2r, &v2v[VLEN * (2 * (sp + 4 * 1))]);
          load_vec(pg, yt2i, &v2v[VLEN * (1 + 2 * (sp + 4 * 1))]);
          load_vec(pg, yt3r, &v2v[VLEN * (2 * (sp + 4 * 2))]);
          load_vec(pg, yt3i, &v2v[VLEN * (1 + 2 * (sp + 4 * 2))]);
          add_vec(pg, yt1r, wt1r);
          add_vec(pg, yt1i, wt1i);
          add_vec(pg, yt2r, wt2r);
          add_vec(pg, yt2i, wt2i);
          add_vec(pg, yt3r, wt3r);
          add_vec(pg, yt3i, wt3i);

          svreal_t xt1r, xt1i, xt2r, xt2i, xt3r, xt3i;
          load_vec(pg, xt1r, &v1[VLEN * (2 * (sp + 4 * 0))]);
          load_vec(pg, xt1i, &v1[VLEN * (1 + 2 * (sp + 4 * 0))]);
          load_vec(pg, xt2r, &v1[VLEN * (2 * (sp + 4 * 1))]);
          load_vec(pg, xt2i, &v1[VLEN * (1 + 2 * (sp + 4 * 1))]);
          load_vec(pg, xt3r, &v1[VLEN * (2 * (sp + 4 * 2))]);
          load_vec(pg, xt3i, &v1[VLEN * (1 + 2 * (sp + 4 * 2))]);

          aypx_vec(pg, a, yt1r, xt1r);
          aypx_vec(pg, a, yt1i, xt1i);
          aypx_vec(pg, a, yt2r, xt2r);
          aypx_vec(pg, a, yt2i, xt2i);
          aypx_vec(pg, a, yt3r, xt3r);
          aypx_vec(pg, a, yt3i, xt3i);

          save_vec(pg, &v2[VLEN * (2 * (sp + 4 * 0))], yt1r);
          save_vec(pg, &v2[VLEN * (1 + 2 * (sp + 4 * 0))], yt1i);
          save_vec(pg, &v2[VLEN * (2 * (sp + 4 * 1))], yt2r);
          save_vec(pg, &v2[VLEN * (1 + 2 * (sp + 4 * 1))], yt2i);
          save_vec(pg, &v2[VLEN * (2 * (sp + 4 * 2))], yt3r);
          save_vec(pg, &v2[VLEN * (1 + 2 * (sp + 4 * 2))], yt3i);
        }
        { // store to sp2, sp4
          svreal_t yt4r, yt4i, yt5r, yt5i, yt6r, yt6i;
          load_vec(pg, yt4r, &v2v[VLEN * (2 * (sp + 1 + 4 * 0))]);
          load_vec(pg, yt4i, &v2v[VLEN * (1 + 2 * (sp + 1 + 4 * 0))]);
          load_vec(pg, yt5r, &v2v[VLEN * (2 * (sp + 1 + 4 * 1))]);
          load_vec(pg, yt5i, &v2v[VLEN * (1 + 2 * (sp + 1 + 4 * 1))]);
          load_vec(pg, yt6r, &v2v[VLEN * (2 * (sp + 1 + 4 * 2))]);
          load_vec(pg, yt6i, &v2v[VLEN * (1 + 2 * (sp + 1 + 4 * 2))]);
          add_vec(pg, yt4r, wt4r);
          add_vec(pg, yt4i, wt4i);
          add_vec(pg, yt5r, wt5r);
          add_vec(pg, yt5i, wt5i);
          add_vec(pg, yt6r, wt6r);
          add_vec(pg, yt6i, wt6i);

          svreal_t xt4r, xt4i, xt5r, xt5i, xt6r, xt6i;
          load_vec(pg, xt4r, &v1[VLEN * (2 * (sp + 1 + 4 * 0))]);
          load_vec(pg, xt4i, &v1[VLEN * (1 + 2 * (sp + 1 + 4 * 0))]);
          load_vec(pg, xt5r, &v1[VLEN * (2 * (sp + 1 + 4 * 1))]);
          load_vec(pg, xt5i, &v1[VLEN * (1 + 2 * (sp + 1 + 4 * 1))]);
          load_vec(pg, xt6r, &v1[VLEN * (2 * (sp + 1 + 4 * 2))]);
          load_vec(pg, xt6i, &v1[VLEN * (1 + 2 * (sp + 1 + 4 * 2))]);

          aypx_vec(pg, a, yt4r, xt4r);
          aypx_vec(pg, a, yt4i, xt4i);
          aypx_vec(pg, a, yt5r, xt5r);
          aypx_vec(pg, a, yt5i, xt5i);
          aypx_vec(pg, a, yt6r, xt6r);
          aypx_vec(pg, a, yt6i, xt6i);

          save_vec(pg, &v2[VLEN * (2 * (sp + 1 + 4 * 0))], yt4r);
          save_vec(pg, &v2[VLEN * (1 + 2 * (sp + 1 + 4 * 0))], yt4i);
          save_vec(pg, &v2[VLEN * (2 * (sp + 1 + 4 * 1))], yt5r);
          save_vec(pg, &v2[VLEN * (1 + 2 * (sp + 1 + 4 * 1))], yt5i);
          save_vec(pg, &v2[VLEN * (2 * (sp + 1 + 4 * 2))], yt6r);
          save_vec(pg, &v2[VLEN * (1 + 2 * (sp + 1 + 4 * 2))], yt6i);
        }
      }
    }
  }


//====================================================================
  template<typename REALTYPE>
  inline void mult_cswinv_chrot(REALTYPE *__restrict v2,
                                REALTYPE *__restrict ct,
                                REALTYPE *__restrict v1)
  {
    svbool_t pg = set_predicate();

    svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i;
    svreal_t vt4r, vt4i, vt5r, vt5i, vt6r, vt6i;
    {
      //  out(0,ic1) = -Uup (0,ic1, id2, ic2) * in2up(id2, ic2);
      //  out(1,ic1) = -Uup (1,ic1, id2, ic2) * in2up(id2, ic2);
      //  out(2,ic1) =  Uup (0,ic1, id2, ic2) * in2up(id2, ic2);
      //  out(3,ic1) =  Uup (1,ic1, id2, ic2) * in2up(id2, ic2);

      set_2sp_up(pg, vt1r, vt1i, v1, 0, 0);
      set_2sp_up(pg, vt2r, vt2i, v1, 0, 1);
      set_2sp_up(pg, vt3r, vt3i, v1, 0, 2);
      set_2sp_up(pg, vt4r, vt4i, v1, 1, 0);
      set_2sp_up(pg, vt5r, vt5i, v1, 1, 1);
      set_2sp_up(pg, vt6r, vt6i, v1, 1, 2);

      svreal_t u11r, u22r, u33r, u44r, u55r, u66r;
      load_vec(pg, u11r, &ct[VLEN * index_clv_up(0, 0, 0, 0, 0)]);
      load_vec(pg, u22r, &ct[VLEN * index_clv_up(0, 1, 0, 1, 0)]);
      load_vec(pg, u33r, &ct[VLEN * index_clv_up(0, 2, 0, 2, 0)]);
      load_vec(pg, u44r, &ct[VLEN * index_clv_up(1, 0, 1, 0, 0)]);
      load_vec(pg, u55r, &ct[VLEN * index_clv_up(1, 1, 1, 1, 0)]);
      load_vec(pg, u66r, &ct[VLEN * index_clv_up(1, 2, 1, 2, 0)]);

      svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i;
      svreal_t wt4r, wt4i, wt5r, wt5i, wt6r, wt6i;
      mul_vec(pg, wt1r, u11r, vt1r); //  = U11 x v1
      mul_vec(pg, wt1i, u11r, vt1i);
      mul_vec(pg, wt2r, u22r, vt2r); // += U22 * v2
      mul_vec(pg, wt2i, u22r, vt2i);
      mul_vec(pg, wt3r, u33r, vt3r); // += U33 x v3
      mul_vec(pg, wt3i, u33r, vt3i);

      svreal_t u12r, u12i;
      load_vec(pg, u12r, &ct[VLEN * index_clv_up(0, 0, 0, 1, 0)]);
      load_vec(pg, u12i, &ct[VLEN * index_clv_up(0, 0, 0, 1, 1)]);

      mul_vec(pg, wt4r, u44r, vt4r); // += U44 x v4
      mul_vec(pg, wt4i, u44r, vt4i);
      mul_vec(pg, wt5r, u55r, vt5r); // += U55 x v5
      mul_vec(pg, wt5i, u55r, vt5i);
      mul_vec(pg, wt6r, u66r, vt6r); // += U66 x v6
      mul_vec(pg, wt6i, u66r, vt6i);

      svreal_t u34r, u34i, u56r, u56i;
      load_vec(pg, u34r, &ct[VLEN * index_clv_up(0, 2, 1, 0, 0)]);
      load_vec(pg, u34i, &ct[VLEN * index_clv_up(0, 2, 1, 0, 1)]);
      load_vec(pg, u56r, &ct[VLEN * index_clv_up(1, 1, 1, 2, 0)]);
      load_vec(pg, u56i, &ct[VLEN * index_clv_up(1, 1, 1, 2, 1)]);

      axpy_vec(pg, wt1r, u12r, vt2r); // += U12 x v2
      axpy_vec(pg, wt1i, u12r, vt2i);
      ymax_vec(pg, wt1r, u12i, vt2i);
      axpy_vec(pg, wt1i, u12i, vt2r);
      axpy_vec(pg, wt2r, u12r, vt1r); // = U21 x v1 = (U12)^* x v1
      axpy_vec(pg, wt2i, u12r, vt1i);
      axpy_vec(pg, wt2r, u12i, vt1i);
      ymax_vec(pg, wt2i, u12i, vt1r);

      svreal_t u23r, u23i;
      load_vec(pg, u23r, &ct[VLEN * index_clv_up(0, 1, 0, 2, 0)]);
      load_vec(pg, u23i, &ct[VLEN * index_clv_up(0, 1, 0, 2, 1)]);
      axpy_vec(pg, wt3r, u34r, vt4r); // += U34 x v4
      axpy_vec(pg, wt3i, u34r, vt4i);
      ymax_vec(pg, wt3r, u34i, vt4i);
      axpy_vec(pg, wt3i, u34i, vt4r);
      axpy_vec(pg, wt4r, u34r, vt3r); // += U43 x v3 = (U34)^* x v3
      axpy_vec(pg, wt4i, u34r, vt3i);
      axpy_vec(pg, wt4r, u34i, vt3i);
      ymax_vec(pg, wt4i, u34i, vt3r);

      svreal_t u45r, u45i;
      load_vec(pg, u45r, &ct[VLEN * index_clv_up(1, 0, 1, 1, 0)]);
      load_vec(pg, u45i, &ct[VLEN * index_clv_up(1, 0, 1, 1, 1)]);
      axpy_vec(pg, wt5r, u56r, vt6r); // += U56 x v6
      axpy_vec(pg, wt5i, u56r, vt6i);
      ymax_vec(pg, wt5r, u56i, vt6i);
      axpy_vec(pg, wt5i, u56i, vt6r);
      axpy_vec(pg, wt6r, u56r, vt5r); // += U65 x v5 = (U56)^* x v5
      axpy_vec(pg, wt6i, u56r, vt5i);
      axpy_vec(pg, wt6r, u56i, vt5i);
      ymax_vec(pg, wt6i, u56i, vt5r);

      svreal_t u13r, u13i;
      load_vec(pg, u13r, &ct[VLEN * index_clv_up(0, 0, 0, 2, 0)]);
      load_vec(pg, u13i, &ct[VLEN * index_clv_up(0, 0, 0, 2, 1)]);
      axpy_vec(pg, wt2r, u23r, vt3r); // += U23 x v3
      axpy_vec(pg, wt2i, u23r, vt3i);
      ymax_vec(pg, wt2r, u23i, vt3i);
      axpy_vec(pg, wt2i, u23i, vt3r);
      axpy_vec(pg, wt3r, u23r, vt2r); // += U32 x v2 = (U23)^* x v2
      axpy_vec(pg, wt3i, u23r, vt2i);
      axpy_vec(pg, wt3r, u23i, vt2i);
      ymax_vec(pg, wt3i, u23i, vt2r);

      svreal_t u24r, u24i;
      load_vec(pg, u24r, &ct[VLEN * index_clv_up(0, 1, 1, 0, 0)]);
      load_vec(pg, u24i, &ct[VLEN * index_clv_up(0, 1, 1, 0, 1)]);
      axpy_vec(pg, wt4r, u45r, vt5r); // += U45 x v5
      axpy_vec(pg, wt4i, u45r, vt5i);
      ymax_vec(pg, wt4r, u45i, vt5i);
      axpy_vec(pg, wt4i, u45i, vt5r);
      axpy_vec(pg, wt5r, u45r, vt4r); // += U54 x v4 = (U45)^* x v4
      axpy_vec(pg, wt5i, u45r, vt4i);
      axpy_vec(pg, wt5r, u45i, vt4i);
      ymax_vec(pg, wt5i, u45i, vt4r);

      svreal_t u35r, u35i;
      load_vec(pg, u35r, &ct[VLEN * index_clv_up(0, 2, 1, 1, 0)]);
      load_vec(pg, u35i, &ct[VLEN * index_clv_up(0, 2, 1, 1, 1)]);
      axpy_vec(pg, wt1r, u13r, vt3r); // += U13 x v3
      axpy_vec(pg, wt1i, u13r, vt3i);
      ymax_vec(pg, wt1r, u13i, vt3i);
      axpy_vec(pg, wt1i, u13i, vt3r);
      axpy_vec(pg, wt3r, u13r, vt1r); // = U31 x v1 = (U13)^* x v1
      axpy_vec(pg, wt3i, u13r, vt1i);
      axpy_vec(pg, wt3r, u13i, vt1i);
      ymax_vec(pg, wt3i, u13i, vt1r);

      svreal_t u46r, u46i;
      load_vec(pg, u46r, &ct[VLEN * index_clv_up(1, 0, 1, 2, 0)]);
      load_vec(pg, u46i, &ct[VLEN * index_clv_up(1, 0, 1, 2, 1)]);
      axpy_vec(pg, wt2r, u24r, vt4r); // += U24 x v4
      axpy_vec(pg, wt2i, u24r, vt4i);
      ymax_vec(pg, wt2r, u24i, vt4i);
      axpy_vec(pg, wt2i, u24i, vt4r);
      axpy_vec(pg, wt4r, u24r, vt2r); // += U42 x v2 = (U24)^* x v2
      axpy_vec(pg, wt4i, u24r, vt2i);
      axpy_vec(pg, wt4r, u24i, vt2i);
      ymax_vec(pg, wt4i, u24i, vt2r);

      svreal_t u15r, u15i;
      load_vec(pg, u15r, &ct[VLEN * index_clv_up(0, 0, 1, 1, 0)]);
      load_vec(pg, u15i, &ct[VLEN * index_clv_up(0, 0, 1, 1, 1)]);
      axpy_vec(pg, wt3r, u35r, vt5r); // += U35 x v5
      axpy_vec(pg, wt3i, u35r, vt5i);
      ymax_vec(pg, wt3r, u35i, vt5i);
      axpy_vec(pg, wt3i, u35i, vt5r);
      axpy_vec(pg, wt5r, u35r, vt3r); // += U53 x v3 = (U35)^* x v3
      axpy_vec(pg, wt5i, u35r, vt3i);
      axpy_vec(pg, wt5r, u35i, vt3i);
      ymax_vec(pg, wt5i, u35i, vt3r);


      svreal_t u26r, u26i;
      load_vec(pg, u26r, &ct[VLEN * index_clv_up(0, 1, 1, 2, 0)]);
      load_vec(pg, u26i, &ct[VLEN * index_clv_up(0, 1, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u46r, vt4r); // += U64 x v4 = (U46)^* x v4
      axpy_vec(pg, wt6i, u46r, vt4i);
      axpy_vec(pg, wt6r, u46i, vt4i);
      ymax_vec(pg, wt6i, u46i, vt4r);
      axpy_vec(pg, wt4r, u46r, vt6r); // += U46 x v6
      axpy_vec(pg, wt4i, u46r, vt6i);
      ymax_vec(pg, wt4r, u46i, vt6i);
      axpy_vec(pg, wt4i, u46i, vt6r);

      svreal_t u14r, u14i;
      load_vec(pg, u14r, &ct[VLEN * index_clv_up(0, 0, 1, 0, 0)]);
      load_vec(pg, u14i, &ct[VLEN * index_clv_up(0, 0, 1, 0, 1)]);
      axpy_vec(pg, wt1r, u15r, vt5r); // += U15 x v5
      axpy_vec(pg, wt1i, u15r, vt5i);
      ymax_vec(pg, wt1r, u15i, vt5i);
      axpy_vec(pg, wt1i, u15i, vt5r);
      axpy_vec(pg, wt5r, u15r, vt1r); // = U51 x v1 = (U15)^* x v1
      axpy_vec(pg, wt5i, u15r, vt1i);
      axpy_vec(pg, wt5r, u15i, vt1i);
      ymax_vec(pg, wt5i, u15i, vt1r);

      svreal_t u36r, u36i;
      load_vec(pg, u36r, &ct[VLEN * index_clv_up(0, 2, 1, 2, 0)]);
      load_vec(pg, u36i, &ct[VLEN * index_clv_up(0, 2, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u26r, vt2r); // += U62 x v2 = (U26)^* x v2
      axpy_vec(pg, wt6i, u26r, vt2i);
      axpy_vec(pg, wt6r, u26i, vt2i);
      ymax_vec(pg, wt6i, u26i, vt2r);
      axpy_vec(pg, wt2r, u26r, vt6r); // += U26 x v6
      axpy_vec(pg, wt2i, u26r, vt6i);
      ymax_vec(pg, wt2r, u26i, vt6i);
      axpy_vec(pg, wt2i, u26i, vt6r);

      svreal_t u25r, u25i;
      load_vec(pg, u25r, &ct[VLEN * index_clv_up(0, 1, 1, 1, 0)]);
      load_vec(pg, u25i, &ct[VLEN * index_clv_up(0, 1, 1, 1, 1)]);
      axpy_vec(pg, wt4r, u14r, vt1r); // = U41 x v1 = (U14)^* x v1
      axpy_vec(pg, wt4i, u14r, vt1i);
      axpy_vec(pg, wt4r, u14i, vt1i);
      ymax_vec(pg, wt4i, u14i, vt1r);
      axpy_vec(pg, wt1r, u14r, vt4r); // += U14 x v4
      axpy_vec(pg, wt1i, u14r, vt4i);
      ymax_vec(pg, wt1r, u14i, vt4i);
      axpy_vec(pg, wt1i, u14i, vt4r);

      svreal_t u16r, u16i;
      load_vec(pg, u16r, &ct[VLEN * index_clv_up(0, 0, 1, 2, 0)]);
      load_vec(pg, u16i, &ct[VLEN * index_clv_up(0, 0, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u36r, vt3r); // += U63 x v3 = (U36)^* x v3
      axpy_vec(pg, wt6i, u36r, vt3i);
      axpy_vec(pg, wt6r, u36i, vt3i);
      ymax_vec(pg, wt6i, u36i, vt3r);
      axpy_vec(pg, wt3r, u36r, vt6r); // += U36 x v6
      axpy_vec(pg, wt3i, u36r, vt6i);
      ymax_vec(pg, wt3r, u36i, vt6i);
      axpy_vec(pg, wt3i, u36i, vt6r);


      axpy_vec(pg, wt2r, u25r, vt5r); // += U25 x v5
      axpy_vec(pg, wt2i, u25r, vt5i);
      ymax_vec(pg, wt2r, u25i, vt5i);
      axpy_vec(pg, wt2i, u25i, vt5r);
      axpy_vec(pg, wt5r, u25r, vt2r); // += U52 x v2 = (U25)^* x v2
      axpy_vec(pg, wt5i, u25r, vt2i);
      axpy_vec(pg, wt5r, u25i, vt2i);
      ymax_vec(pg, wt5i, u25i, vt2r);

      axpy_vec(pg, wt1r, u16r, vt6r); // += U16 x v6
      axpy_vec(pg, wt1i, u16r, vt6i);
      ymax_vec(pg, wt1r, u16i, vt6i);
      axpy_vec(pg, wt1i, u16i, vt6r);
      axpy_vec(pg, wt6r, u16r, vt1r); // = U61 x v1 = (U16)^* x v1
      axpy_vec(pg, wt6i, u16r, vt1i);
      axpy_vec(pg, wt6r, u16i, vt1i);
      ymax_vec(pg, wt6i, u16i, vt1r);

      save_vec(pg, &v2[VLEN * (2 * (2 + 4 * 0))], wt1r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (2 + 4 * 0))], wt1i);
      save_vec(pg, &v2[VLEN * (2 * (2 + 4 * 1))], wt2r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (2 + 4 * 1))], wt2i);
      save_vec(pg, &v2[VLEN * (2 * (2 + 4 * 2))], wt3r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (2 + 4 * 2))], wt3i);

      save_vec(pg, &v2[VLEN * (2 * (3 + 4 * 0))], wt4r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (3 + 4 * 0))], wt4i);
      save_vec(pg, &v2[VLEN * (2 * (3 + 4 * 1))], wt5r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (3 + 4 * 1))], wt5i);
      save_vec(pg, &v2[VLEN * (2 * (3 + 4 * 2))], wt6r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (3 + 4 * 2))], wt6i);

      flip_sign(pg, wt1r);
      flip_sign(pg, wt1i);
      flip_sign(pg, wt2r);
      flip_sign(pg, wt2i);
      flip_sign(pg, wt3r);
      flip_sign(pg, wt3i);
      flip_sign(pg, wt4r);
      flip_sign(pg, wt4i);
      flip_sign(pg, wt5r);
      flip_sign(pg, wt5i);
      flip_sign(pg, wt6r);
      flip_sign(pg, wt6i);

      save_vec(pg, &v2[VLEN * (2 * (0 + 4 * 0))], wt1r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (0 + 4 * 0))], wt1i);
      save_vec(pg, &v2[VLEN * (2 * (0 + 4 * 1))], wt2r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (0 + 4 * 1))], wt2i);
      save_vec(pg, &v2[VLEN * (2 * (0 + 4 * 2))], wt3r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (0 + 4 * 2))], wt3i);

      save_vec(pg, &v2[VLEN * (2 * (1 + 4 * 0))], wt4r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (1 + 4 * 0))], wt4i);
      save_vec(pg, &v2[VLEN * (2 * (1 + 4 * 1))], wt5r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (1 + 4 * 1))], wt5i);
      save_vec(pg, &v2[VLEN * (2 * (1 + 4 * 2))], wt6r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (1 + 4 * 2))], wt6i);
    }


    {
      set_2sp_dn(pg, vt1r, vt1i, v1, 0, 0);
      set_2sp_dn(pg, vt2r, vt2i, v1, 0, 1);
      set_2sp_dn(pg, vt3r, vt3i, v1, 0, 2);
      set_2sp_dn(pg, vt4r, vt4i, v1, 1, 0);
      set_2sp_dn(pg, vt5r, vt5i, v1, 1, 1);
      set_2sp_dn(pg, vt6r, vt6i, v1, 1, 2);

      svreal_t u11r, u22r, u33r, u44r, u55r, u66r;
      load_vec(pg, u11r, &ct[VLEN * index_clv_dn(0, 0, 0, 0, 0)]);
      load_vec(pg, u22r, &ct[VLEN * index_clv_dn(0, 1, 0, 1, 0)]);
      load_vec(pg, u33r, &ct[VLEN * index_clv_dn(0, 2, 0, 2, 0)]);
      load_vec(pg, u44r, &ct[VLEN * index_clv_dn(1, 0, 1, 0, 0)]);
      load_vec(pg, u55r, &ct[VLEN * index_clv_dn(1, 1, 1, 1, 0)]);
      load_vec(pg, u66r, &ct[VLEN * index_clv_dn(1, 2, 1, 2, 0)]);

      svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i;
      svreal_t wt4r, wt4i, wt5r, wt5i, wt6r, wt6i;
      mul_vec(pg, wt1r, u11r, vt1r); //  = U11 x v1
      mul_vec(pg, wt1i, u11r, vt1i);
      mul_vec(pg, wt2r, u22r, vt2r); // += U22 * v2
      mul_vec(pg, wt2i, u22r, vt2i);
      mul_vec(pg, wt3r, u33r, vt3r); // += U33 x v3
      mul_vec(pg, wt3i, u33r, vt3i);

      svreal_t u12r, u12i;
      load_vec(pg, u12r, &ct[VLEN * index_clv_dn(0, 0, 0, 1, 0)]);
      load_vec(pg, u12i, &ct[VLEN * index_clv_dn(0, 0, 0, 1, 1)]);

      mul_vec(pg, wt4r, u44r, vt4r); // += U44 x v4
      mul_vec(pg, wt4i, u44r, vt4i);
      mul_vec(pg, wt5r, u55r, vt5r); // += U55 x v5
      mul_vec(pg, wt5i, u55r, vt5i);
      mul_vec(pg, wt6r, u66r, vt6r); // += U66 x v6
      mul_vec(pg, wt6i, u66r, vt6i);

      svreal_t u34r, u34i, u56r, u56i;
      load_vec(pg, u34r, &ct[VLEN * index_clv_dn(0, 2, 1, 0, 0)]);
      load_vec(pg, u34i, &ct[VLEN * index_clv_dn(0, 2, 1, 0, 1)]);
      load_vec(pg, u56r, &ct[VLEN * index_clv_dn(1, 1, 1, 2, 0)]);
      load_vec(pg, u56i, &ct[VLEN * index_clv_dn(1, 1, 1, 2, 1)]);

      axpy_vec(pg, wt1r, u12r, vt2r); // += U12 x v2
      axpy_vec(pg, wt1i, u12r, vt2i);
      ymax_vec(pg, wt1r, u12i, vt2i);
      axpy_vec(pg, wt1i, u12i, vt2r);
      axpy_vec(pg, wt2r, u12r, vt1r); // = U21 x v1 = (U12)^* x v1
      axpy_vec(pg, wt2i, u12r, vt1i);
      axpy_vec(pg, wt2r, u12i, vt1i);
      ymax_vec(pg, wt2i, u12i, vt1r);

      svreal_t u23r, u23i;
      load_vec(pg, u23r, &ct[VLEN * index_clv_dn(0, 1, 0, 2, 0)]);
      load_vec(pg, u23i, &ct[VLEN * index_clv_dn(0, 1, 0, 2, 1)]);
      axpy_vec(pg, wt3r, u34r, vt4r); // += U34 x v4
      axpy_vec(pg, wt3i, u34r, vt4i);
      ymax_vec(pg, wt3r, u34i, vt4i);
      axpy_vec(pg, wt3i, u34i, vt4r);
      axpy_vec(pg, wt4r, u34r, vt3r); // += U43 x v3 = (U34)^* x v3
      axpy_vec(pg, wt4i, u34r, vt3i);
      axpy_vec(pg, wt4r, u34i, vt3i);
      ymax_vec(pg, wt4i, u34i, vt3r);

      svreal_t u45r, u45i;
      load_vec(pg, u45r, &ct[VLEN * index_clv_dn(1, 0, 1, 1, 0)]);
      load_vec(pg, u45i, &ct[VLEN * index_clv_dn(1, 0, 1, 1, 1)]);
      axpy_vec(pg, wt5r, u56r, vt6r); // += U56 x v6
      axpy_vec(pg, wt5i, u56r, vt6i);
      ymax_vec(pg, wt5r, u56i, vt6i);
      axpy_vec(pg, wt5i, u56i, vt6r);
      axpy_vec(pg, wt6r, u56r, vt5r); // += U65 x v5 = (U56)^* x v5
      axpy_vec(pg, wt6i, u56r, vt5i);
      axpy_vec(pg, wt6r, u56i, vt5i);
      ymax_vec(pg, wt6i, u56i, vt5r);

      svreal_t u13r, u13i;
      load_vec(pg, u13r, &ct[VLEN * index_clv_dn(0, 0, 0, 2, 0)]);
      load_vec(pg, u13i, &ct[VLEN * index_clv_dn(0, 0, 0, 2, 1)]);
      axpy_vec(pg, wt2r, u23r, vt3r); // += U23 x v3
      axpy_vec(pg, wt2i, u23r, vt3i);
      ymax_vec(pg, wt2r, u23i, vt3i);
      axpy_vec(pg, wt2i, u23i, vt3r);
      axpy_vec(pg, wt3r, u23r, vt2r); // += U32 x v2 = (U23)^* x v2
      axpy_vec(pg, wt3i, u23r, vt2i);
      axpy_vec(pg, wt3r, u23i, vt2i);
      ymax_vec(pg, wt3i, u23i, vt2r);

      svreal_t u24r, u24i;
      load_vec(pg, u24r, &ct[VLEN * index_clv_dn(0, 1, 1, 0, 0)]);
      load_vec(pg, u24i, &ct[VLEN * index_clv_dn(0, 1, 1, 0, 1)]);
      axpy_vec(pg, wt4r, u45r, vt5r); // += U45 x v5
      axpy_vec(pg, wt4i, u45r, vt5i);
      ymax_vec(pg, wt4r, u45i, vt5i);
      axpy_vec(pg, wt4i, u45i, vt5r);
      axpy_vec(pg, wt5r, u45r, vt4r); // += U54 x v4 = (U45)^* x v4
      axpy_vec(pg, wt5i, u45r, vt4i);
      axpy_vec(pg, wt5r, u45i, vt4i);
      ymax_vec(pg, wt5i, u45i, vt4r);

      svreal_t u35r, u35i;
      load_vec(pg, u35r, &ct[VLEN * index_clv_dn(0, 2, 1, 1, 0)]);
      load_vec(pg, u35i, &ct[VLEN * index_clv_dn(0, 2, 1, 1, 1)]);
      axpy_vec(pg, wt1r, u13r, vt3r); // += U13 x v3
      axpy_vec(pg, wt1i, u13r, vt3i);
      ymax_vec(pg, wt1r, u13i, vt3i);
      axpy_vec(pg, wt1i, u13i, vt3r);
      axpy_vec(pg, wt3r, u13r, vt1r); // = U31 x v1 = (U13)^* x v1
      axpy_vec(pg, wt3i, u13r, vt1i);
      axpy_vec(pg, wt3r, u13i, vt1i);
      ymax_vec(pg, wt3i, u13i, vt1r);

      svreal_t u46r, u46i;
      load_vec(pg, u46r, &ct[VLEN * index_clv_dn(1, 0, 1, 2, 0)]);
      load_vec(pg, u46i, &ct[VLEN * index_clv_dn(1, 0, 1, 2, 1)]);
      axpy_vec(pg, wt2r, u24r, vt4r); // += U24 x v4
      axpy_vec(pg, wt2i, u24r, vt4i);
      ymax_vec(pg, wt2r, u24i, vt4i);
      axpy_vec(pg, wt2i, u24i, vt4r);
      axpy_vec(pg, wt4r, u24r, vt2r); // += U42 x v2 = (U24)^* x v2
      axpy_vec(pg, wt4i, u24r, vt2i);
      axpy_vec(pg, wt4r, u24i, vt2i);
      ymax_vec(pg, wt4i, u24i, vt2r);

      svreal_t u15r, u15i;
      load_vec(pg, u15r, &ct[VLEN * index_clv_dn(0, 0, 1, 1, 0)]);
      load_vec(pg, u15i, &ct[VLEN * index_clv_dn(0, 0, 1, 1, 1)]);
      axpy_vec(pg, wt3r, u35r, vt5r); // += U35 x v5
      axpy_vec(pg, wt3i, u35r, vt5i);
      ymax_vec(pg, wt3r, u35i, vt5i);
      axpy_vec(pg, wt3i, u35i, vt5r);
      axpy_vec(pg, wt5r, u35r, vt3r); // += U53 x v3 = (U35)^* x v3
      axpy_vec(pg, wt5i, u35r, vt3i);
      axpy_vec(pg, wt5r, u35i, vt3i);
      ymax_vec(pg, wt5i, u35i, vt3r);

      svreal_t u26r, u26i;
      load_vec(pg, u26r, &ct[VLEN * index_clv_dn(0, 1, 1, 2, 0)]);
      load_vec(pg, u26i, &ct[VLEN * index_clv_dn(0, 1, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u46r, vt4r); // += U64 x v4 = (U46)^* x v4
      axpy_vec(pg, wt6i, u46r, vt4i);
      axpy_vec(pg, wt6r, u46i, vt4i);
      ymax_vec(pg, wt6i, u46i, vt4r);
      axpy_vec(pg, wt4r, u46r, vt6r); // += U46 x v6
      axpy_vec(pg, wt4i, u46r, vt6i);
      ymax_vec(pg, wt4r, u46i, vt6i);
      axpy_vec(pg, wt4i, u46i, vt6r);

      svreal_t u14r, u14i;
      load_vec(pg, u14r, &ct[VLEN * index_clv_dn(0, 0, 1, 0, 0)]);
      load_vec(pg, u14i, &ct[VLEN * index_clv_dn(0, 0, 1, 0, 1)]);
      axpy_vec(pg, wt1r, u15r, vt5r); // += U15 x v5
      axpy_vec(pg, wt1i, u15r, vt5i);
      ymax_vec(pg, wt1r, u15i, vt5i);
      axpy_vec(pg, wt1i, u15i, vt5r);
      axpy_vec(pg, wt5r, u15r, vt1r); // = U51 x v1 = (U15)^* x v1
      axpy_vec(pg, wt5i, u15r, vt1i);
      axpy_vec(pg, wt5r, u15i, vt1i);
      ymax_vec(pg, wt5i, u15i, vt1r);

      svreal_t u36r, u36i;
      load_vec(pg, u36r, &ct[VLEN * index_clv_dn(0, 2, 1, 2, 0)]);
      load_vec(pg, u36i, &ct[VLEN * index_clv_dn(0, 2, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u26r, vt2r); // += U62 x v2 = (U26)^* x v2
      axpy_vec(pg, wt6i, u26r, vt2i);
      axpy_vec(pg, wt6r, u26i, vt2i);
      ymax_vec(pg, wt6i, u26i, vt2r);
      axpy_vec(pg, wt2r, u26r, vt6r); // += U26 x v6
      axpy_vec(pg, wt2i, u26r, vt6i);
      ymax_vec(pg, wt2r, u26i, vt6i);
      axpy_vec(pg, wt2i, u26i, vt6r);

      svreal_t u25r, u25i;
      load_vec(pg, u25r, &ct[VLEN * index_clv_dn(0, 1, 1, 1, 0)]);
      load_vec(pg, u25i, &ct[VLEN * index_clv_dn(0, 1, 1, 1, 1)]);
      axpy_vec(pg, wt4r, u14r, vt1r); // = U41 x v1 = (U14)^* x v1
      axpy_vec(pg, wt4i, u14r, vt1i);
      axpy_vec(pg, wt4r, u14i, vt1i);
      ymax_vec(pg, wt4i, u14i, vt1r);
      axpy_vec(pg, wt1r, u14r, vt4r); // += U14 x v4
      axpy_vec(pg, wt1i, u14r, vt4i);
      ymax_vec(pg, wt1r, u14i, vt4i);
      axpy_vec(pg, wt1i, u14i, vt4r);

      svreal_t u16r, u16i;
      load_vec(pg, u16r, &ct[VLEN * index_clv_dn(0, 0, 1, 2, 0)]);
      load_vec(pg, u16i, &ct[VLEN * index_clv_dn(0, 0, 1, 2, 1)]);
      axpy_vec(pg, wt6r, u36r, vt3r); // += U63 x v3 = (U36)^* x v3
      axpy_vec(pg, wt6i, u36r, vt3i);
      axpy_vec(pg, wt6r, u36i, vt3i);
      ymax_vec(pg, wt6i, u36i, vt3r);
      axpy_vec(pg, wt3r, u36r, vt6r); // += U36 x v6
      axpy_vec(pg, wt3i, u36r, vt6i);
      ymax_vec(pg, wt3r, u36i, vt6i);
      axpy_vec(pg, wt3i, u36i, vt6r);


      axpy_vec(pg, wt2r, u25r, vt5r); // += U25 x v5
      axpy_vec(pg, wt2i, u25r, vt5i);
      ymax_vec(pg, wt2r, u25i, vt5i);
      axpy_vec(pg, wt2i, u25i, vt5r);
      axpy_vec(pg, wt5r, u25r, vt2r); // += U52 x v2 = (U25)^* x v2
      axpy_vec(pg, wt5i, u25r, vt2i);
      axpy_vec(pg, wt5r, u25i, vt2i);
      ymax_vec(pg, wt5i, u25i, vt2r);

      axpy_vec(pg, wt1r, u16r, vt6r); // += U16 x v6
      axpy_vec(pg, wt1i, u16r, vt6i);
      ymax_vec(pg, wt1r, u16i, vt6i);
      axpy_vec(pg, wt1i, u16i, vt6r);
      axpy_vec(pg, wt6r, u16r, vt1r); // = U61 x v1 = (U16)^* x v1
      axpy_vec(pg, wt6i, u16r, vt1i);
      axpy_vec(pg, wt6r, u16i, vt1i);
      ymax_vec(pg, wt6i, u16i, vt1r);

      svreal_t xt1r, xt1i, xt2r, xt2i, xt3r, xt3i;
      svreal_t xt4r, xt4i, xt5r, xt5i, xt6r, xt6i;


      load_vec(pg, xt1r, &v2[VLEN * (2 * (0 + 4 * 0))]);
      load_vec(pg, xt1i, &v2[VLEN * (1 + 2 * (0 + 4 * 0))]);
      load_vec(pg, xt2r, &v2[VLEN * (2 * (0 + 4 * 1))]);
      load_vec(pg, xt2i, &v2[VLEN * (1 + 2 * (0 + 4 * 1))]);
      load_vec(pg, xt3r, &v2[VLEN * (2 * (0 + 4 * 2))]);
      load_vec(pg, xt3i, &v2[VLEN * (1 + 2 * (0 + 4 * 2))]);
      load_vec(pg, xt4r, &v2[VLEN * (2 * (1 + 4 * 0))]);
      load_vec(pg, xt4i, &v2[VLEN * (1 + 2 * (1 + 4 * 0))]);
      load_vec(pg, xt5r, &v2[VLEN * (2 * (1 + 4 * 1))]);
      load_vec(pg, xt5i, &v2[VLEN * (1 + 2 * (1 + 4 * 1))]);
      load_vec(pg, xt6r, &v2[VLEN * (2 * (1 + 4 * 2))]);
      load_vec(pg, xt6i, &v2[VLEN * (1 + 2 * (1 + 4 * 2))]);

      add_vec(pg, xt1r, wt1r);
      add_vec(pg, xt1i, wt1i);
      add_vec(pg, xt2r, wt2r);
      add_vec(pg, xt2i, wt2i);
      add_vec(pg, xt3r, wt3r);
      add_vec(pg, xt3i, wt3i);
      add_vec(pg, xt4r, wt4r);
      add_vec(pg, xt4i, wt4i);
      add_vec(pg, xt5r, wt5r);
      add_vec(pg, xt5i, wt5i);
      add_vec(pg, xt6r, wt6r);
      add_vec(pg, xt6i, wt6i);

      save_vec(pg, &v2[VLEN * (2 * (0 + 4 * 0))], xt1r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (0 + 4 * 0))], xt1i);
      save_vec(pg, &v2[VLEN * (2 * (0 + 4 * 1))], xt2r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (0 + 4 * 1))], xt2i);
      save_vec(pg, &v2[VLEN * (2 * (0 + 4 * 2))], xt3r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (0 + 4 * 2))], xt3i);
      save_vec(pg, &v2[VLEN * (2 * (1 + 4 * 0))], xt4r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (1 + 4 * 0))], xt4i);
      save_vec(pg, &v2[VLEN * (2 * (1 + 4 * 1))], xt5r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (1 + 4 * 1))], xt5i);
      save_vec(pg, &v2[VLEN * (2 * (1 + 4 * 2))], xt6r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (1 + 4 * 2))], xt6i);

      load_vec(pg, xt1r, &v2[VLEN * (2 * (2 + 4 * 0))]);
      load_vec(pg, xt1i, &v2[VLEN * (1 + 2 * (2 + 4 * 0))]);
      load_vec(pg, xt2r, &v2[VLEN * (2 * (2 + 4 * 1))]);
      load_vec(pg, xt2i, &v2[VLEN * (1 + 2 * (2 + 4 * 1))]);
      load_vec(pg, xt3r, &v2[VLEN * (2 * (2 + 4 * 2))]);
      load_vec(pg, xt3i, &v2[VLEN * (1 + 2 * (2 + 4 * 2))]);
      load_vec(pg, xt4r, &v2[VLEN * (2 * (3 + 4 * 0))]);
      load_vec(pg, xt4i, &v2[VLEN * (1 + 2 * (3 + 4 * 0))]);
      load_vec(pg, xt5r, &v2[VLEN * (2 * (3 + 4 * 1))]);
      load_vec(pg, xt5i, &v2[VLEN * (1 + 2 * (3 + 4 * 1))]);
      load_vec(pg, xt6r, &v2[VLEN * (2 * (3 + 4 * 2))]);
      load_vec(pg, xt6i, &v2[VLEN * (1 + 2 * (3 + 4 * 2))]);

      add_vec(pg, xt1r, wt1r);
      add_vec(pg, xt1i, wt1i);
      add_vec(pg, xt2r, wt2r);
      add_vec(pg, xt2i, wt2i);
      add_vec(pg, xt3r, wt3r);
      add_vec(pg, xt3i, wt3i);
      add_vec(pg, xt4r, wt4r);
      add_vec(pg, xt4i, wt4i);
      add_vec(pg, xt5r, wt5r);
      add_vec(pg, xt5i, wt5i);
      add_vec(pg, xt6r, wt6r);
      add_vec(pg, xt6i, wt6i);

      save_vec(pg, &v2[VLEN * (2 * (2 + 4 * 0))], xt1r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (2 + 4 * 0))], xt1i);
      save_vec(pg, &v2[VLEN * (2 * (2 + 4 * 1))], xt2r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (2 + 4 * 1))], xt2i);
      save_vec(pg, &v2[VLEN * (2 * (2 + 4 * 2))], xt3r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (2 + 4 * 2))], xt3i);
      save_vec(pg, &v2[VLEN * (2 * (3 + 4 * 0))], xt4r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (3 + 4 * 0))], xt4i);
      save_vec(pg, &v2[VLEN * (2 * (3 + 4 * 1))], xt5r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (3 + 4 * 1))], xt5i);
      save_vec(pg, &v2[VLEN * (2 * (3 + 4 * 2))], xt6r);
      save_vec(pg, &v2[VLEN * (1 + 2 * (3 + 4 * 2))], xt6i);
    }
  }


//====================================================================
} // nameless namespace end

#endif
