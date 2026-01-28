/*!
        @file    aprojection_Stout_SU3-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: aoyama $
        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$
        @version $LastChangedRevision: 2314 $
*/

#include "lib/Smear/projection_Stout_SU3.h"

// The following implementation only valid for Nc = 3 case.
#define  NC     3
#define  NVC    6
#define  NDF   18


template<typename AFIELD>
const std::string AProjection_Stout_SU3<AFIELD>::class_name
                                    = "AProjection_Stout_SU3<AFIELD>";

//====================================================================
namespace {

  inline int idxr(int ic1, int ic2){ return    2*(ic1 + ic2 * NC); }
  inline int idxi(int ic1, int ic2){ return 1+ 2*(ic1 + ic2 * NC); }

  inline void acos_vec(svbool_t& pg, svreal_t& x, const svreal_t& y){
    Vsimd_t xt, yt;
    save_vec(pg, &yt, y);
    for(int k = 0; k < VLEN; ++k){
      xt.v[k] = acos(yt.v[k]);
    }
    load_vec(pg, x, &xt.v[0]);
  }

  inline void cos_vec(svbool_t& pg, svreal_t& x, const svreal_t& y){
    Vsimd_t xt, yt;
    save_vec(pg, &yt, y);
    for(int k = 0; k < VLEN; ++k){
      xt.v[k] = cos(yt.v[k]);
    }
    load_vec(pg, x, &xt.v[0]);
  }

  inline void sin_vec(svbool_t& pg, svreal_t& x, const svreal_t& y){
    Vsimd_t xt, yt;
    save_vec(pg, &yt, y);
    for(int k = 0; k < VLEN; ++k){
      xt.v[k] = sin(yt.v[k]);
    }
    load_vec(pg, x, &xt.v[0]);
  }

  inline void set_unit(Vsimd_t* u){
    clear_vec(u, NDF);
    set_vec(&u[ 0], real_t(1.0), 1);
    set_vec(&u[ 8], real_t(1.0), 1);
    set_vec(&u[16], real_t(1.0), 1);
  }

  //  inline void load_u(Vsimd_t* u, real_t* up){
  //    for(int idf = 0; idf < NDF; ++idf){
  //      load_vec(&u[idf], &up[VLEN * idf], NDF);
  //    }
  //  }

  inline void mult_nn(Vsimd_t* u1, Vsimd_t* u2, Vsimd_t* u3){
    svbool_t pg = set_predicate();
    for(int ic1 = 0; ic1 < NC; ++ic1){
      for(int ic2 = 0; ic2 < NC; ++ic2){
        svreal_t utr, uti;
        set_vec(pg, utr, real_t(0.0));
        set_vec(pg, uti, real_t(0.0));
        for(int k = 0; k < NC; ++k){
          svreal_t u2r, u2i, u3r, u3i;
          load_vec(pg, u2r, &u2[idxr(ic1,k)].v[0]);
          load_vec(pg, u2i, &u2[idxi(ic1,k)].v[0]);
          load_vec(pg, u3r, &u3[idxr(k,ic2)].v[0]);
          load_vec(pg, u3i, &u3[idxi(k,ic2)].v[0]);
          axpy_vec(pg, utr, u2r, u3r);
          ymax_vec(pg, utr, u2i, u3i);
          axpy_vec(pg, uti, u2r, u3i);
          axpy_vec(pg, uti, u2i, u3r);
	}
        save_vec(pg, &u1[idxr(ic1,ic2)], utr);
        save_vec(pg, &u1[idxi(ic1,ic2)], uti);
      }
    }
  }

  inline void mult_nd(Vsimd_t* u1, Vsimd_t* u2, Vsimd_t* u3){
    svbool_t pg = set_predicate();
    for(int ic1 = 0; ic1 < NC; ++ic1){
      for(int ic2 = 0; ic2 < NC; ++ic2){
        svreal_t utr, uti;
        set_vec(pg, utr, real_t(0.0));
        set_vec(pg, uti, real_t(0.0));
        for(int k = 0; k < NC; ++k){
          svreal_t u2r, u2i, u3r, u3i;
          load_vec(pg, u2r, &u2[idxr(ic1,k)].v[0]);
          load_vec(pg, u2i, &u2[idxi(ic1,k)].v[0]);
          load_vec(pg, u3r, &u3[idxr(ic2,k)].v[0]);
          load_vec(pg, u3i, &u3[idxi(ic2,k)].v[0]);
	  axpy_vec(pg, utr, u2r, u3r);
          axpy_vec(pg, utr, u2i, u3i);
          ymax_vec(pg, uti, u2r, u3i);
          axpy_vec(pg, uti, u2i, u3r);
	}
        save_vec(pg, &u1[idxr(ic1,ic2)], utr);
        save_vec(pg, &u1[idxi(ic1,ic2)], uti);
      }
    }
  }

  inline void mat_at(Vsimd_t* u){
    svbool_t pg = set_predicate();
    svreal_t tri;
    set_vec(pg, tri, real_t(0.0));
    for(int ic1 = 0; ic1 < NC; ++ic1){
      for(int ic2 = ic1; ic2 < NC; ++ic2){
        if(ic1 != ic2){
          svreal_t ur, ui, unr, uni, utr, uti;	
          load_vec(pg, unr, &u[idxr(ic1,ic2)].v[0]);
	  load_vec(pg, uni, &u[idxi(ic1,ic2)].v[0]);
	  load_vec(pg, utr, &u[idxr(ic2,ic1)].v[0]);
	  load_vec(pg, uti, &u[idxi(ic2,ic1)].v[0]);
	  sub_vec(pg, ur, unr, utr);
          add_vec(pg, ui, uni, uti);
          scal_vec(pg, ur, real_t(0.5));
          scal_vec(pg, ui, real_t(0.5));
	  save_vec(pg, &u[idxr(ic1,ic2)], ur);
          save_vec(pg, &u[idxi(ic1,ic2)], ui);
          scal_vec(pg, ur, real_t(-1.0));
	  save_vec(pg, &u[idxr(ic2,ic1)], ur);
          save_vec(pg, &u[idxi(ic2,ic1)], ui);
	}else{
          svreal_t uti;
	  load_vec(pg, uti, &u[idxi(ic1,ic1)].v[0]);
          add_vec(pg, tri, uti);
	}
      }
    }
    svreal_t utr, uti;
    set_vec(pg, utr, 0.0);
    for(int ic1 = 0; ic1 < NC; ++ic1){
      load_vec(pg, uti, &u[idxi(ic1,ic1)].v[0]);
      axpy_vec(pg, uti, -1.0/3.0, tri);
      save_vec(pg, &u[idxr(ic1,ic1)], utr);
      save_vec(pg, &u[idxi(ic1,ic1)], uti);
    }
  }

  inline void func_xi0(svreal_t& v, const svreal_t& w)
  {
    Vsimd_t vt, wt;
    svbool_t pg = set_predicate();
    save_vec(pg, &wt, w);
    for(int k = 0; k < VLEN; ++k){
      if (wt.v[k] == 0.0) {
        vt.v[k] = 1.0;
      } else {
        vt.v[k] = sin(wt.v[k]) / wt.v[k];
      }
    }
    load_vec(pg, v, &vt.v[0]);
  }

  inline void set_fj2(Vsimd_t *f0, Vsimd_t *f1, Vsimd_t *f2,
                      const Vsimd_t& u, const Vsimd_t& w)
  {
    // f0, f1, f2 are complex vectors

    // const double xi0   = func_xi0(w);
    // const double u2    = u * u;
    // const double w2    = w * w;
    // const double cos_w = cos(w);
    svbool_t pg = set_predicate();
    svreal_t xi0, u2, w2, cos_w, ut, wt;
    load_vec(pg, ut, &u.v[0]);
    load_vec(pg, wt, &w.v[0]);
    func_xi0(xi0, wt);
    mul_vec(pg, u2, ut, ut);
    mul_vec(pg, w2, wt, wt);
    cos_vec(pg, cos_w, wt);

    // const double cos_u = cos(u);
    // const double sin_u = sin(u);
    svreal_t cos_u, sin_u;
    cos_vec(pg, cos_u, ut);
    sin_vec(pg, sin_u, ut);

    // const dcomplex emiu = cmplx(cos_u, -sin_u);
    // const dcomplex e2iu = cmplx(cos_u * cos_u - sin_u * sin_u,
    //                             2.0 * sin_u * cos_u);
    svreal_t emiu_r, emiu_i, e2iu_r, e2iu_i;
    set_vec(pg, emiu_r, cos_u);
    set_vec(pg, emiu_i, real_t(-1.0), sin_u);
    mul_vec(pg, e2iu_r, cos_u, cos_u);
    ymax_vec(pg, e2iu_r, sin_u, sin_u);
    mul_vec(pg, e2iu_i, sin_u, cos_u);
    scal_vec(pg, e2iu_i, real_t(2.0));

    //const dcomplex h0 =  e2iu * cmplx(u2 - w2, 0.0)
    //                   + emiu * cmplx(8.0 * u2 * cos_w,
    //                                 2.0 * u * (3.0 * u2 + w2) * xi0);
    svreal_t h0_r, h0_i, vt1, vt2, vt3;
    sub_vec(pg, vt1, u2, w2);
    mul_vec(pg, h0_r, e2iu_r, vt1);
    mul_vec(pg, h0_i, e2iu_i, vt1);
    mul_vec(pg, vt1, u2, cos_w);
    scal_vec(pg, vt1, real_t(8.0));
    axpy_vec(pg, h0_r, emiu_r, vt1);
    axpy_vec(pg, h0_i, emiu_i, vt1);
    set_vec(pg, vt1, 3.0, u2);
    add_vec(pg, vt2, vt1, w2);
    mul_vec(pg, vt3, ut, xi0);
    mul_vec(pg, vt1, vt2, vt3);
    mul_vec(pg, vt2, emiu_r, vt1);
    mul_vec(pg, vt3, emiu_i, vt1);
    axpy_vec(pg, h0_r, real_t(-2.0), vt3);
    axpy_vec(pg, h0_i, real_t( 2.0), vt2);

    //const dcomplex h1 =  e2iu * cmplx(2.0 * u, 0.0)
    //                   - emiu * cmplx(2.0 * u * cos_w,
    //                                  -(3.0 * u2 - w2) * xi0);
    svreal_t h1_r, h1_i;
    mul_vec(pg, h1_r, e2iu_r, ut);
    mul_vec(pg, h1_i, e2iu_i, ut);
    scal_vec(pg, h1_r, 2.0);
    scal_vec(pg, h1_i, 2.0);
    mul_vec(pg, vt1, ut, cos_w);
    mul_vec(pg, vt2, emiu_r, vt1);
    mul_vec(pg, vt3, emiu_i, vt1);
    axpy_vec(pg, h1_r, -2.0, vt2);
    axpy_vec(pg, h1_i, -2.0, vt3);
    set_vec(pg, vt1, 3.0, u2);
    sub_vec(pg, vt2, vt1, w2);
    mul_vec(pg, vt3, vt2, xi0);
    mul_vec(pg, vt1, emiu_r, vt3);
    mul_vec(pg, vt2, emiu_i, vt3);
    axpy_vec(pg, h1_r, real_t(-1.0), vt2);
    axpy_vec(pg, h1_i, real_t( 1.0), vt1);

    // const dcomplex h2 = e2iu - emiu * cmplx(cos_w, 3.0 * u * xi0);
    svreal_t h2_r, h2_i;
    set_vec(pg, h2_r, e2iu_r);
    set_vec(pg, h2_i, e2iu_i);
    ymax_vec(pg, h2_r, emiu_r, cos_w);
    ymax_vec(pg, h2_i, emiu_i, cos_w);
    mul_vec(pg, vt3, ut, xi0);
    mul_vec(pg, vt1, emiu_r, vt3);
    mul_vec(pg, vt2, emiu_i, vt3);
    axpy_vec(pg, h2_r, real_t( 3.0), vt2);
    axpy_vec(pg, h2_i, real_t(-3.0), vt1);

    //  const double fden = 1.0 / (9.0 * u2 - w2);
    set_vec(pg, vt1, real_t(9.0), u2);
    sub_vec(pg, vt2, vt1, w2);
    set_vec(pg, vt3, 1.0);
    div_vec(pg, vt1, vt3, vt2);

    // f0 = h0 * fden;
    // f1 = h1 * fden;
    // f2 = h2 * fden;
    svreal_t ft_r, ft_i;
    mul_vec(pg, ft_r, h0_r, vt1);
    mul_vec(pg, ft_i, h0_i, vt1);
    save_vec(pg, &f0[0], ft_r);
    save_vec(pg, &f0[1], ft_i);

    mul_vec(pg, ft_r, h1_r, vt1);
    mul_vec(pg, ft_i, h1_i, vt1);
    save_vec(pg, &f1[0], ft_r);
    save_vec(pg, &f1[1], ft_i);

    mul_vec(pg, ft_r, h2_r, vt1);
    mul_vec(pg, ft_i, h2_i, vt1);
    save_vec(pg, &f2[0], ft_r);
    save_vec(pg, &f2[1], ft_i);
  }

  inline void set_uw2(Vsimd_t& u, Vsimd_t& w,
                      const Vsimd_t* iQ2, const Vsimd_t* iQ3)
  {
    //  real_t c0 = -(iQ3[idxi(0,0)] + iQ3[idxi(1,1)] + iQ3[idxi(2, 2)]) / 3.0;
    //  real_t c1 = -0.5 * (iQ2[idxr(0, 0)] + iQ2[idxr(1, 1)] + iQ2[idxr(2, 2)]);
    //  real_t c13r  = sqrt(c1/3.0);
    //  real_t c0max = 2.0 * c13r * c13r * c13r;

    svbool_t pg = set_predicate();

    svreal_t c0, c1, c13r, c0max, theta, ct, ct2;

    svreal_t iq3_00, iq3_11, iq3_22;
    load_vec(pg, iq3_00, &iQ3[idxi(0,0)].v[0]);
    load_vec(pg, iq3_11, &iQ3[idxi(1,1)].v[0]);
    load_vec(pg, iq3_22, &iQ3[idxi(2,2)].v[0]);

    add_vec(pg, ct, iq3_00, iq3_11);
    add_vec(pg, c0, ct, iq3_22);
    scal_vec(pg, c0, real_t(-1.0/3.0));
 
    svreal_t iq2_00, iq2_11, iq2_22;
    load_vec(pg, iq2_00, &iQ2[idxr(0,0)].v[0]);
    load_vec(pg, iq2_11, &iQ2[idxr(1,1)].v[0]);
    load_vec(pg, iq2_22, &iQ2[idxr(2,2)].v[0]);

    add_vec(pg, ct, iq2_00, iq2_11);
    add_vec(pg, c1, ct, iq2_22);
    scal_vec(pg, c1, real_t(-0.5));

    set_vec(pg, ct, real_t(1.0/3.0), c1);
    sqrt_vec(pg, c13r, ct);

    mul_vec(pg, ct, c13r, c13r);
    mul_vec(pg, c0max, ct, c13r);
    scal_vec(pg, c0max, real_t(2.0));

    //  real_t theta = acos(c0/c0max);
    div_vec(pg, ct, c0, c0max);
    acos_vec(pg, theta, ct);
  
    //- output
    //  u = c13r * cos(theta / 3.0);
    //  w = sqrt(c1) * sin(theta / 3.0);
    scal_vec(pg, theta, real_t(1.0/3.0));
    cos_vec(pg, ct, theta);
    mul_vec(pg, ct2, c13r, ct);
    save_vec(pg, &u, ct2);

    sin_vec(pg, ct, theta);
    sqrt_vec(pg, ct2, c1);
    svreal_t ct3;
    mul_vec(pg, ct3, ct2, ct);
    save_vec(pg, &w, ct3);
  }

  inline void set_expiQ(Vsimd_t *e_iQ,
                        const Vsimd_t *f0,  const Vsimd_t *f1,
			const Vsimd_t *f2,  const Vsimd_t *iQ0,
			const Vsimd_t *iQ1,  const Vsimd_t *iQ2)
  {
    svbool_t pg = set_predicate();
    svreal_t f0r, f0i, f1r, f1i, f2r, f2i;
    load_vec(pg, f0r, &f0[0].v[0]);
    load_vec(pg, f0i, &f0[1].v[0]);
    load_vec(pg, f1r, &f1[0].v[0]);
    load_vec(pg, f1i, &f1[1].v[0]);
    load_vec(pg, f2r, &f2[0].v[0]);
    load_vec(pg, f2i, &f2[1].v[0]);

    for(int ic1 = 0; ic1 < NC; ++ic1){
      for(int ic2 = 0; ic2 < NC; ++ic2){
        // dcomplex qt =   f0 * cmplx(iQ0.r(cc),  iQ0.i(cc))
        //               + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
        //               - f2 * cmplx(iQ2.r(cc),  iQ2.i(cc));
        // e_iQ.set_r(cc, real(qt));
        // e_iQ.set_i(cc, imag(qt));
        svreal_t iq0r, iq0i, iq1r, iq1i, iq2r, iq2i;
        load_vec(pg, iq0r, &iQ0[idxr(ic1,ic2)].v[0]);
        load_vec(pg, iq0i, &iQ0[idxi(ic1,ic2)].v[0]);
        load_vec(pg, iq1r, &iQ1[idxr(ic1,ic2)].v[0]);
        load_vec(pg, iq1i, &iQ1[idxi(ic1,ic2)].v[0]);
        load_vec(pg, iq2r, &iQ2[idxr(ic1,ic2)].v[0]);
        load_vec(pg, iq2i, &iQ2[idxi(ic1,ic2)].v[0]);

        svreal_t ex_r, ex_i;
        mul_vec( pg, ex_r, f0r, iq0r);
        ymax_vec(pg, ex_r, f0i, iq0i);
        axpy_vec(pg, ex_r, f1r, iq1i);
        axpy_vec(pg, ex_r, f1i, iq1r);
        ymax_vec(pg, ex_r, f2r, iq2r);
        axpy_vec(pg, ex_r, f2i, iq2i);
        save_vec(pg, &e_iQ[idxr(ic1,ic2)], ex_r);
	  
        mul_vec( pg, ex_i, f0r, iq0i);
        axpy_vec(pg, ex_i, f0i, iq0r);
        ymax_vec(pg, ex_i, f1r, iq1r);
        axpy_vec(pg, ex_i, f1i, iq1i);
        ymax_vec(pg, ex_i, f2r, iq2i);
        ymax_vec(pg, ex_i, f2i, iq2r);
        save_vec(pg, &e_iQ[idxi(ic1,ic2)], ex_i);
      }
    }
  }

}

//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::init()
{
  assert(CommonParameters::Nc() == NC);

  m_vl = CommonParameters::Vlevel();

  // strict check
  if (CommonParameters::Nc() != NC) {
    vout.crucial(m_vl, "Error at %s: Nc = 3 is needed, but Nc = %d\n",
                 class_name.c_str(), CommonParameters::Nc());
    exit(EXIT_FAILURE);
  }

  m_Ndf = 2 * NC * NC;
  m_Nst  = CommonParameters::Nvol();
  m_Nstv = m_Nst/VLEN;
  int Ndim = CommonParameters::Ndim();

  m_v1.reset(m_Ndf, m_Nst, Ndim);
  m_v2.reset(m_Ndf, m_Nst, Ndim);
  m_v3.reset(m_Ndf, m_Nst, Ndim);

  m_flop = 0;
  m_time = 0.0;

}


//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  return;
}


//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::print_stat()
{
  const double gflops = 1.0e-9 * double(m_flop) / m_time;

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  total time: %f\n", m_time);
  vout.general(m_vl, "  total flop: %d\n", m_flop);
  vout.general(m_vl, "  GFlops    : %f\n", gflops);
}


//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::project(Field_G& U,
                                            const double alpha,
                                            const Field_G& Cst,
                                            const Field_G& Uorg)
{
#pragma omp barrier

  AIndex_lex<real_t,AFIELD::IMPL> index;

  convert_gauge(index, m_v1, Cst);
  convert_gauge(index, m_v2, Uorg);

  project(m_v3, alpha, m_v1, m_v2);

  reverse_gauge(index, U, m_v3);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::project(AFIELD& U,
                                            const double alpha,
                                            const AFIELD& Cst,
                                            const AFIELD& Uorg)
{
  // in stout projection, parameter alpha is dummy.

#pragma omp barrier

  const double time0 = Communicator::get_time();

  const int Nex  = Uorg.nex();
  const int Nvol = Uorg.nvol();
  const int NinG = Uorg.nin();

  assert(Cst.nex() == Nex);
  assert(Cst.nvol() == Nvol);
  assert(U.nex() == Nex);
  assert(U.nvol() == Nvol);

  real_t *cst  = const_cast<AFIELD*>(&Cst)->ptr(0);
  real_t *uorg = const_cast<AFIELD*>(&Uorg)->ptr(0);
  real_t *up   = U.ptr(0);

  Vsimd_t iQ0[NDF];
  set_unit(iQ0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

  AIndex_lex<real_t,QXS> index;

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = is; site < ns; ++site) {

      Vsimd_t ut[NDF];
      load_vec(ut, &uorg[VLEN*NDF*(site + m_Nstv * mu)], NDF);
 
      Vsimd_t ct[NDF];
      load_vec(ct, &cst[VLEN*NDF*(site + m_Nstv * mu)], NDF);

      Vsimd_t iQ1[NDF];
      mult_nd(iQ1, ct, ut);
      mat_at(iQ1);

      Vsimd_t iQ2[NDF];
      mult_nn(iQ2, iQ1, iQ1);

      Vsimd_t iQ3[NDF];
      mult_nn(iQ3, iQ1, iQ2);

      Vsimd_t e_iQ[NDF];

      real_t norm;
      norm2_vec(norm, iQ1, NDF);
      // this evaluates summation of |iQ1|^2 over sites in a vector.
      // better definition is norm2 on site by site.

      if (norm > 1.0e-10) {
        Vsimd_t u, w;
        set_uw2(u, w, iQ2, iQ3);
        Vsimd_t f0[2], f1[2], f2[2];
        set_fj2(f0, f1, f2, u, w);
	set_expiQ(e_iQ, f0, f1, f2, iQ0, iQ1, iQ2);
      } else {
        set_unit(e_iQ);
      }

      Vsimd_t ut2[NDF];
      mult_nn(ut2, e_iQ, ut);

      save_vec(&up[VLEN*NDF*(site + m_Nstv * mu)], ut2, NDF);
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::project_alt(Field_G& U,
                                            const double alpha,
                                            const Field_G& Cst,
                                            const Field_G& Uorg)
{
#pragma omp barrier

  const double time0 = Communicator::get_time();

  // in stout projection, parameter alpha is dummy.

  const int Nex  = Uorg.nex();
  const int Nvol = Uorg.nvol();
  const int NinG = Uorg.nin();

  assert(Cst.nex() == Nex);
  assert(Cst.nvol() == Nvol);
  assert(U.nex() == Nex);
  assert(U.nvol() == Nvol);

  Mat_SU_N iQ0(NC);
  iQ0.unit();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = is; site < ns; ++site) {
      Mat_SU_N ut(NC);
      Uorg.mat(ut, site, mu);

      Mat_SU_N ct(NC);
      Cst.mat(ct, site, mu);

      Mat_SU_N iQ1(NC);
      iQ1.mult_nd(ct, ut);
      iQ1.at();

      Mat_SU_N iQ2(NC);
      iQ2.mult_nn(iQ1, iQ1);

      Mat_SU_N iQ3(NC);
      iQ3.mult_nn(iQ1, iQ2);

      Mat_SU_N e_iQ(NC);

      double norm = iQ1.norm2();
      if (norm > 1.0e-10) {
        double u, w;
        set_uw(u, w, iQ2, iQ3);

        dcomplex f0, f1, f2;
        set_fj(f0, f1, f2, u, w);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt =   f0 * cmplx(iQ0.r(cc),  iQ0.i(cc))
                        + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - f2 * cmplx(iQ2.r(cc),  iQ2.i(cc));
          e_iQ.set_r(cc, real(qt));
          e_iQ.set_i(cc, imag(qt));
        }
      } else {
        //  vout.general(m_vl,"project: |iQ1|^2 too small: %lf. Set e_iQ=1.\n",norm);
        e_iQ.unit();
      }

      Mat_SU_N ut2(NC);
      ut2.mult_nn(e_iQ, ut);
      U.set_mat(site, mu, ut2);
    }
  }

  /*
  unsigned long count;
  double time;
  KEK_FopCountFinish(id,&count,&time);
  m_time += time;
  m_flop += count;
  */
  const double time1 = Communicator::get_time();
  m_time += time1 - time0;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::exp_iQ(Field_G& e_iQ, const Field_G& iQ)
{
#pragma omp barrier

  const int Nvol = iQ.nvol();
  const int Nex  = iQ.nex();

  Mat_SU_N iQ0(NC);

  iQ0.unit();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = is; site < ns; ++site) {
      Mat_SU_N iQ1 = iQ.mat(site, mu);
      Mat_SU_N iQ2 = iQ1 * iQ1;
      Mat_SU_N iQ3 = iQ1 * iQ2;

      double norm = iQ1.norm2();
      if (norm > 1.0e-10) {
        double u, w;
        set_uw(u, w, iQ2, iQ3);

        dcomplex f0, f1, f2;
        set_fj(f0, f1, f2, u, w);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = f0 * cmplx(iQ0.r(cc), iQ0.i(cc))
                        + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - f2 * cmplx(iQ2.r(cc), iQ2.i(cc));
          e_iQ.set_ri(cc, site, mu, real(qt), imag(qt));
        }
      } else {
        //      vout.general(m_vl,"exp_iQ: |iQ1|^2 too small: %lf. Set e_iQ=1.\n",norm);
        e_iQ.set_mat(site, mu, iQ0);
      }
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::force_recursive(Field_G& Xi,
                                           Field_G& iTheta,
                                           const double alpha,
                                           const Field_G& Sigmap,
                                           const Field_G& Cst,
                                           const Field_G& Uorg)
{
#pragma omp barrier
  // in stout projection, parameter alpha is dummy.

  //  int id = 31;
  //  KEK_FopCountStart(id);
  const double time0 = Communicator::get_time();

  const int Nvol = CommonParameters::Nvol();
  const int Nex = Xi.nex();

  assert(Xi.nvol() == Nvol);
  assert(iTheta.nvol() == Nvol);
  assert(Sigmap.nvol() == Nvol);
  assert(Cst.nvol() == Nvol);
  assert(Uorg.nvol() == Nvol);
  assert(iTheta.nex() == Nex);
  assert(Sigmap.nex() == Nex);
  assert(Cst.nex() == Nex);
  assert(Uorg.nex() == Nex);

  Mat_SU_N iQ0(NC);
  iQ0.unit();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int mu = 0; mu < Nex; ++mu) {
    for (int site = is; site < ns; ++site) {
      //! C_tmp \f$=C_\mu(x)\f$
      Mat_SU_N C_tmp(NC);
      Cst.mat(C_tmp, site, mu);

      //! U_tmp \f$=U_\mu(x)\f$
      Mat_SU_N U_tmp(NC);
      Uorg.mat(U_tmp, site, mu);

      // Sigmap_tmp \f$=\Sigma_\mu'(x)\f$
      Mat_SU_N Sigmap_tmp(NC);
      Sigmap.mat(Sigmap_tmp, site, mu);

      //! iQ1 \f$=iQ_\mu\f$
      Mat_SU_N iQ1(NC);
      iQ1.mult_nd(C_tmp, U_tmp);
      iQ1.at();

      Mat_SU_N iQ2(NC);
      iQ2.mult_nn(iQ1, iQ1);

      Mat_SU_N iQ3(NC);
      iQ3.mult_nn(iQ1, iQ2);

      // In order to aviod 1Q1=0
      Mat_SU_N e_iQ(NC), iGamma(NC);

      double norm = iQ1.norm2();
      if (norm > 1.0e-10) {
        double u, w;
        set_uw(u, w, iQ2, iQ3);

        dcomplex f0, f1, f2;
        set_fj(f0, f1, f2, u, w);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = f0 * cmplx(iQ0.r(cc), iQ0.i(cc))
                        + f1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - f2 * cmplx(iQ2.r(cc), iQ2.i(cc));
          e_iQ.set(cc, real(qt), imag(qt));
        }

        double xi0   = func_xi0(w);
        double xi1   = func_xi1(w);
        double u2    = u * u;
        double w2    = w * w;
        double cos_w = cos(w);

        dcomplex emiu = cmplx(cos(u), -sin(u));
        dcomplex e2iu = cmplx(cos(2.0 * u), sin(2.0 * u));

        dcomplex r01 = cmplx(2.0 * u, 2.0 * (u2 - w2)) * e2iu
                       + emiu * cmplx(16.0 * u * cos_w + 2.0 * u * (3.0 * u2 + w2) * xi0,
                                      -8.0 * u2 * cos_w + 2.0 * (9.0 * u2 + w2) * xi0);

        dcomplex r11 = cmplx(2.0, 4.0 * u) * e2iu
                       + emiu * cmplx(-2.0 * cos_w + (3.0 * u2 - w2) * xi0,
                                      2.0 * u * cos_w + 6.0 * u * xi0);

        dcomplex r21 = cmplx(0.0, 2.0) * e2iu
                       + emiu * cmplx(-3.0 * u * xi0, cos_w - 3.0 * xi0);

        dcomplex r02 = cmplx(-2.0, 0.0) * e2iu
                       + emiu * cmplx(-8.0 * u2 * xi0,
                                      2.0 * u * (cos_w + xi0 + 3.0 * u2 * xi1));

        dcomplex r12 = emiu * cmplx(2.0 * u * xi0,
                                    -cos_w - xi0 + 3.0 * u2 * xi1);

        dcomplex r22 = emiu * cmplx(xi0, -3.0 * u * xi1);

        double fden = 1.0 / (2 * (9.0 * u2 - w2) * (9.0 * u2 - w2));

        dcomplex b10 = cmplx(2.0 * u, 0.0) * r01 + cmplx(3.0 * u2 - w2, 0.0) * r02
                       - cmplx(30.0 * u2 + 2.0 * w2, 0.0) * f0;
        dcomplex b11 = cmplx(2.0 * u, 0.0) * r11 + cmplx(3.0 * u2 - w2, 0.0) * r12
                       - cmplx(30.0 * u2 + 2.0 * w2, 0.0) * f1;
        dcomplex b12 = cmplx(2.0 * u, 0.0) * r21 + cmplx(3.0 * u2 - w2, 0.0) * r22
                       - cmplx(30.0 * u2 + 2.0 * w2, 0.0) * f2;

        dcomplex b20 = r01 - cmplx(3.0 * u, 0.0) * r02 - cmplx(24.0 * u, 0.0) * f0;
        dcomplex b21 = r11 - cmplx(3.0 * u, 0.0) * r12 - cmplx(24.0 * u, 0.0) * f1;
        dcomplex b22 = r21 - cmplx(3.0 * u, 0.0) * r22 - cmplx(24.0 * u, 0.0) * f2;

        b10 *= cmplx(fden, 0.0);
        b11 *= cmplx(fden, 0.0);
        b12 *= cmplx(fden, 0.0);
        b20 *= cmplx(fden, 0.0);
        b21 *= cmplx(fden, 0.0);
        b22 *= cmplx(fden, 0.0);

        Mat_SU_N B1(NC), B2(NC);
        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt1 = b10 * cmplx(iQ0.r(cc), iQ0.i(cc))
                         + b11 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                         - b12 * cmplx(iQ2.r(cc), iQ2.i(cc));
          B1.set(cc, real(qt1), imag(qt1));

          dcomplex qt2 = b20 * cmplx(iQ0.r(cc), iQ0.i(cc))
                         + b21 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                         - b22 * cmplx(iQ2.r(cc), iQ2.i(cc));
          B2.set(cc, real(qt2), imag(qt2));
        }

        Mat_SU_N USigmap(NC);
        USigmap.mult_nn(U_tmp, Sigmap_tmp);

        Mat_SU_N tmp1(NC);
        tmp1.mult_nn(USigmap, B1);

        Mat_SU_N tmp2(NC);
        tmp2.mult_nn(USigmap, B2);

        dcomplex tr1 = cmplx(tmp1.r(0) + tmp1.r(4) + tmp1.r(8),
                             tmp1.i(0) + tmp1.i(4) + tmp1.i(8));
        dcomplex tr2 = cmplx(tmp2.r(0) + tmp2.r(4) + tmp2.r(8),
                             tmp2.i(0) + tmp2.i(4) + tmp2.i(8));

        Mat_SU_N iQUS(NC);
        iQUS.mult_nn(iQ1, USigmap);

        Mat_SU_N iUSQ(NC);
        iUSQ.mult_nn(USigmap, iQ1);

        for (int cc = 0; cc < NC * NC; ++cc) {
          dcomplex qt = tr1 * cmplx(iQ1.i(cc), -iQ1.r(cc))
                        - tr2 * cmplx(iQ2.r(cc), iQ2.i(cc))
                        + f1 * cmplx(USigmap.r(cc), USigmap.i(cc))
                        + f2 * cmplx(iQUS.i(cc), -iQUS.r(cc))
                        + f2 * cmplx(iUSQ.i(cc), -iUSQ.r(cc));
          iGamma.set(cc, -imag(qt), real(qt));
        }
      } else {
        // vout.general(m_vl,"force_recursive: |iQ1|^2 too small: %lf. Set e_iQ=1.\n",norm);
        iGamma.zero();
        e_iQ.unit();
      }

      //! iGamma \f$=i\Lambda\f$
      iGamma.at();

      Mat_SU_N iTheta_tmp(NC);
      iTheta_tmp.mult_nn(iGamma, U_tmp);

      //! iTheta \f$=i\Lambda U_\mu(x)\f$
      iTheta.set_mat(site, mu, iTheta_tmp);

      Mat_SU_N Xi_tmp(NC);
      Xi_tmp.mult_nn(Sigmap_tmp, e_iQ);
      Xi_tmp.multadd_dn(C_tmp, iGamma);

      //! Xi \f$=\Sigma_\mu'(x)\exp(iQ_\mu(x))+C_\mu^\dagger i\Lambda_\mu(x)\f$.
      Xi.set_mat(site, mu, Xi_tmp);
    }
  }

  /*
  unsigned long count;
  double time;
  KEK_FopCountFinish(id,&count,&time);
  m_time += time;
  m_flop += count;
  */
  const double time1 = Communicator::get_time();
  m_time += time1 - time0;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::set_fj(dcomplex& f0, dcomplex& f1, dcomplex& f2,
                                  const double& u, const double& w)
{
  const double xi0   = func_xi0(w);
  const double u2    = u * u;
  const double w2    = w * w;
  const double cos_w = cos(w);

  const double cos_u = cos(u);
  const double sin_u = sin(u);

  const dcomplex emiu = cmplx(cos_u, -sin_u);
  const dcomplex e2iu = cmplx(cos_u * cos_u - sin_u * sin_u,
                              2.0 * sin_u * cos_u);

  const dcomplex h0 =  e2iu * cmplx(u2 - w2, 0.0)
                     + emiu * cmplx(8.0 * u2 * cos_w,
                                   2.0 * u * (3.0 * u2 + w2) * xi0);
  const dcomplex h1 =  e2iu * cmplx(2 * u, 0.0)
                     - emiu * cmplx(2.0 * u * cos_w,
                                    -(3.0 * u2 - w2) * xi0);
  const dcomplex h2 = e2iu - emiu * cmplx(cos_w, 3.0 * u * xi0);

  const double fden = 1.0 / (9.0 * u2 - w2);

  f0 = h0 * fden;
  f1 = h1 * fden;
  f2 = h2 * fden;

}

//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::set_uw(double& u, double& w,
                                  const Mat_SU_N& iQ2, const Mat_SU_N& iQ3)
{
  const double c0    = -(iQ3.i(0, 0) + iQ3.i(1, 1) + iQ3.i(2, 2)) / 3.0;
  const double c1    = -0.5 * (iQ2.r(0, 0) + iQ2.r(1, 1) + iQ2.r(2, 2));
  const double c13r  = sqrt(c1 / 3.0);
  const double c0max = 2.0 * c13r * c13r * c13r;

  const double theta = acos(c0 / c0max);

  //- output
  u = c13r * cos(theta / 3.0);
  w = sqrt(c1) * sin(theta / 3.0);
}


//====================================================================
template<typename AFIELD>
double AProjection_Stout_SU3<AFIELD>::func_xi0(const double w)
{
  if (w == 0.0) {
    return 1.0;
  } else {
    return sin(w) / w;
  }
}


//====================================================================
template<typename AFIELD>
double AProjection_Stout_SU3<AFIELD>::func_xi1(const double w)
{
  if (w < 0.25) {
    const double        w2 = w * w;
    const static double c0 = -1.0 / 3.0;
    const static double c1 =  1.0 / 30.0;
    const static double c2 = -1.0 / 840.0;
    const static double c3 =  1.0 / 45360.0;
    const static double c4 = -1.0 / 3991680.0;

    return c0 + w2 * (c1 + w2 * (c2 + w2 * (c3 + w2 * c4)));
  } else {
    return (w * cos(w) - sin(w)) / (w * w * w);
  }
}


//====================================================================
template<typename AFIELD>
void AProjection_Stout_SU3<AFIELD>::exp_iQ_bf(Field_G& e_iQ, const Field_G& iQ)
{
#pragma omp barrier
  // brute force version of exponentiation: for check

  const static int Nprec = 32;

  const int Nvol = iQ.nvol();
  const int Nex  = iQ.nex();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      Mat_SU_N u0(NC);
      u0.unit();

      Mat_SU_N u1(NC);
      u1.unit();

      Mat_SU_N h1 = iQ.mat(site, ex);

      for (int iprec = 0; iprec < Nprec; ++iprec) {
        double   exf = 1.0 / (Nprec - iprec);
        Mat_SU_N u2  = h1 * u1;

        u2 *= exf;
        u1  = u2;
        u1 += u0;
      }

      u1.reunit();
      e_iQ.set_mat(site, ex, u1);
    }
  }

#pragma omp barrier
}

//============================================================END=====
