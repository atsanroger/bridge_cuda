/*!
        @file    afopr_Overlap_5d.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-05 20:48:53 #$

        @version $LastChangedRevision: 2530 $
*/

#ifndef AFOPR_OVERLAP_5D_INCLUDED
#define AFOPR_OVERLAP_5D_INCLUDED

#include "fopr_Wilson_eo.h"
#include "Tools/math_Sign_Zolotarev.h"
#include "bridge_complex.h"
#include "complexTraits.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! 5-dimensional overlap fermion operator.

/*!
    This class is the 5-dimensional overlap operator used in
    5D overlap solver.
    The implementation is in even-odd site index, and only the
    functionality needed in the 5D solver is ready.
    Kernel fermion operator (even-odd) must be constructed
    outside this class and to be supplied at the construction.
                                     [24 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
    Changed to a template class in ver.2.0.
                                     [12 Feb 2022 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Overlap_5d : public AFopr<AFIELD>
{
 public:
  typedef AFopr<AFIELD>  AFOPR;
  typedef typename AFIELD::real_t real_t;
  typedef typename ComplexTraits<real_t>::complex_t complex_t;

  static const std::string class_name;

 private:

  // input parameters
  real_t m_mq;     //!< quark mass
  real_t m_M0;     //!< domain-wall height
  int m_Np;        //!< number of poles in rational approx.
  real_t m_x_min;  //!< lower range of approximate sign function
  real_t m_x_max;  //!< upper range of approximate sign function
  int m_Niter_ms;         //!< max iteration of shiftsolver (dummy)
  real_t m_Stop_cond_ms;  //!< stopping condition of shift solver (dummy)
  std::vector<int> m_boundary;  //!< boundary consition

  std::string m_kernel_type;  //!< kernel type (if given)
  std::string m_repr;         //!< gamma-matrix repr. (if given)

  Bridge::VerboseLevel m_vl;

  std::string m_mode;

  //- low-mode subtraction
  int m_Nsbt;
  std::vector<real_t> m_ev;
  std::vector<AFIELD> m_vk;

  // local variables
  int m_Nin, m_Nvol2, m_Nex, m_Ndim;

  //- following parameters are used pnly in this class
  std::vector<real_t> m_cl;
  std::vector<real_t> m_bl;
  std::vector<real_t> m_sigma;

  std::vector<real_t> m_p_sqrt;
  std::vector<real_t> m_q_sqrt;
  real_t m_p0_parameter, m_R_parameter, m_h, m_u0;
  std::vector<real_t> m_rl;
  std::vector<real_t> m_sl;

  std::vector<real_t> m_prf;
  std::valarray<complex_t> m_u0c_e, m_u0cinv_e;
  std::valarray<complex_t> m_u0c_o, m_u0cinv_o;

  AFOPR *m_fopr_w;
  bool m_kernel_created;  //!< whether kernel is created in this object
  
  Index_eo m_index_eo;

  AFIELD m_z1, m_z2, m_w1, m_v1, m_v2;
  AFIELD m_t1, m_t2, m_t3;

 public:

  AFopr_Overlap_5d(const Parameters& params) {init(params); }

  AFopr_Overlap_5d(AFOPR *fopr, const Parameters& params)
  {init(fopr, params); }

  DEPRECATED
  AFopr_Overlap_5d(AFOPR *fopr) : m_fopr_w(fopr) { init(); }

  void set_parameters(const Parameters& params);
  void set_parameters(const real_t mq, const real_t M0, const int Np,
                      const real_t x_min, const real_t x_max,
                      const int Niter, const real_t Stop_cond,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U) { m_fopr_w->set_config(U); }

  void set_lowmodes(const int Nsbt, std::vector<real_t> *,
                                    std::vector<AFIELD> *);

  void mult(AFIELD& v, const AFIELD& w);

  void mult_dag(AFIELD& v, const AFIELD& w);

  void mult_gm5(AFIELD& v, const AFIELD& w);

  void set_mode(const std::string mode);

  std::string get_mode() const { return m_mode; }

  void DdagD_eo(AFIELD& v, const AFIELD&);
  void DD_5d_eo(AFIELD& v, const AFIELD& w, const int jd);

  void Mopr_5d_eo(AFIELD& v, const AFIELD& w, const int ieo);

  void LUprecond(AFIELD& v, const AFIELD& w, const int ieo);

  void Proj_H_eo(const int ieo1, const int ieo2,
                 AFIELD& v1, const AFIELD& w1);

  void Proj_L_mult_eo(const int ieo1, const int ieo2,
                      AFIELD& v1, const AFIELD& w1);

  void mult_u0inv(AFIELD& v1, const AFIELD& w1, const int ieo);

  void Calc_Coeff_u0inv();
  void Solv_Coeff_u0inv(std::valarray<complex_t>&,
                        const std::valarray<complex_t>&,
                        const std::valarray<complex_t>&);

  void set_coefficients();

  real_t norm_c(const std::valarray<complex_t>&);

  real_t innerprod_c(const std::valarray<complex_t>&,
                     const std::valarray<complex_t>&);

  void mult_WdagW(std::valarray<complex_t>&,
                  const std::valarray<complex_t>&,
                  const std::valarray<complex_t>&);

  //! returns true if additional field conversion is needed.
  virtual bool needs_convert()
  { return m_fopr_w->needs_convert(); }

  //! converts a Field object into other format if necessary.
  virtual void convert(AFIELD& v, const Field& w);

  //! reverses to a Field object from other format if necessary.
  virtual void reverse(Field& v, const AFIELD& w);

  int field_nin() { return m_fopr_w->field_nin(); }
  int field_nvol(){ return m_fopr_w->field_nvol(); }
  int field_nex() { return 2 * m_Np + 1; }

  //! this returns the number of floating point operations.
  double flop_count();

private:

  void init(const Parameters& params);

  void init(AFOPR *fopr, const Parameters& params);

  void init();

  void tidyup();


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object(AFOPR *fopr)
  {
    return new AFopr_Overlap_5d<AFIELD>(fopr);
  }

  static AFopr<AFIELD> *create_object_with_params(AFOPR *fopr,
                                                  const Parameters& params)
  {
    return new AFopr_Overlap_5d<AFIELD>(fopr, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AFopr<AFIELD>::Factory_fopr::Register("Overlap_5d", create_object);
    init &= AFopr<AFIELD>::Factory_fopr_params::Register("Overlap_5d",
                                                    create_object_with_params);
    return init;
  }
#endif

};
#endif
