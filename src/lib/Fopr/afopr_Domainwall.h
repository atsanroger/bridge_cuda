/*!
      @file    afopr_Domainwall.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2513 $
*/

#ifndef AFOPR_DOMAINWALL_INCLUDED
#define AFOPR_DOMAINWALL_INCLUDED

#include <vector>
#include <string>

#include "lib/Fopr/afopr.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

class Field;
class Field_G;

//! Domain-wall fermion operator.

/*!
    This class implements the domain-wall fermion operator.
    The first version was implemented only for the standard
    (Shamir's) form.                [24 Dec 2011 H.Matsufuru]
    Later it was generalized to cover the Mobius form as a
    template class in the alternative branch.
                                    [26 Feb 2019 H.Matsufuru]
    In ver.2.0, it was merged to trunk by renaming to
    AFopr_Domainwall.
                                    [06 Mar 2022 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Domainwall : public AFopr<AFIELD>
{
 public:
  typedef AFopr<AFIELD>             AFOPR;
  typedef typename AFIELD::real_t   real_t;
  static const std::string class_name;

 private:
  // parameters common to overlap fermion
  real_t m_mq;                 //!< quark mass
  real_t m_M0;                 //!< domain-wall height
  int m_Ns;                    //!< size of fifth-dimension
  std::vector<int> m_boundary; //!< boundary conditions
  std::vector<real_t> m_b;     //!< coefficient b (array)
  std::vector<real_t> m_c;     //!< coefficient c (array)
  Bridge::VerboseLevel m_vl;   //!< verbose level
  std::string m_repr;

  std::string m_mode;

  int m_NinF, m_Nvol, m_Ndim;

  AFOPR *m_foprw;
  std::string m_kernel_type;
  bool m_kernel_created;         //!< whether kernel is created in this object

  AFIELD m_w4, m_v4, m_t4, m_y4; //!< working 4d vectors.
  AFIELD m_w1, m_v1, m_v2;       //!< working 5d vectors.

  // for convert and reverse
  Field m_w4lex;
  AFIELD m_v4lex;

  // for preconditioning
  std::vector<real_t> m_dp;
  std::vector<real_t> m_dm;
  std::vector<real_t> m_e;
  std::vector<real_t> m_f;
  real_t m_g;

 public:
  //! constructor.
  AFopr_Domainwall(AFOPR *fopr) { init(fopr); }

  AFopr_Domainwall(const Parameters& params) { init(params); }

  //! constructor.
  AFopr_Domainwall(AFOPR *fopr, const Parameters& params)
  { init(fopr, params); }

  //! destructor.
  ~AFopr_Domainwall() { tidyup(); }

  void set_parameters(const Parameters& params);

  //! set parameters in the case of Moebius domain-wall.
  void set_parameters(const real_t mq, const real_t M0,
                      const int Ns, const std::vector<int> bc,
                      const real_t b, const real_t c);

  //! set parameters in the case of Moebius domain-wall.
  void set_parameters(const real_t mq, const real_t M0,
                      const int Ns, const std::vector<int> bc,
                      const std::vector<real_t> vec_b,
                      const std::vector<real_t> vec_c);

  void get_parameters(Parameters& params) const;

  //! set parameters of kernel operaotr.
  void set_kernel_parameters(const Parameters& params);

  //! set parameters for preconditioning.
  void set_precond_parameters();

  //! set coefficients if they depend in s.
  void set_coefficients(const std::vector<real_t> b,
                        const std::vector<real_t> c);

  //! this class needs convert of fermion field.
  bool needs_convert() { return m_foprw->needs_convert(); }

  //! convert Field to AField for this class.
  void convert(AFIELD&, const Field&);

  //! reverse AField to Field.
  void reverse(Field&, const AFIELD&);

  void set_config(Field *U) { m_foprw->set_config(U); }

  void set_config(unique_ptr<Field_G>& U)
  { m_foprw->set_config(U.get()); }

  void set_mode(std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(AFIELD& v, const AFIELD& w);

  void mult_dag(AFIELD& v, const AFIELD& w);

  //! mult with specified mode.
  //! Possible modes are: "D", "Ddag", "DdagD", "DDdag",
  //!   "Prec", "D_prec", "Ddag_prec", "DdagD_prec".
  void mult(AFIELD& v, const AFIELD& w, std::string mode);

  //! mult_dag with specified mode.
  void mult_dag(AFIELD& v, const AFIELD& w, std::string mode);

  void mult_gm5(AFIELD&, const AFIELD&);

  void mult_chproj_4d(AFIELD&, const AFIELD&, const int ipm);

  int field_nin() { return m_NinF; }
  int field_nvol() { return m_Nvol; }
  int field_nex() { return m_Ns; }

  //! this returns the number of floating point number operations.
  double flop_count() { return flop_count(m_mode); }

  //! flop-count for specified mode.
  double flop_count(std::string mode);

  void DdagD(AFIELD&, const AFIELD&);
  void DDdag(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);
  void Hdag(AFIELD&, const AFIELD&);
  void mult_gm5R(AFIELD&, const AFIELD&);
  void mult_R(AFIELD&, const AFIELD&);

  //  preconditioner
  void DdagD_prec(AFIELD&, const AFIELD&);
  void D_prec(AFIELD&, const AFIELD&);
  void Ddag_prec(AFIELD&, const AFIELD&);
  void Prec(AFIELD&, const AFIELD&);
  void Precdag(AFIELD&, const AFIELD&);

 private:

  //! initial setup.
  void init(const Parameters& params);

  void init(AFOPR *fopr, const Parameters& params);

  void init(AFOPR *fopr);

  //! final tidyup.
  void tidyup();


  void L_inv(AFIELD&, const AFIELD&);
  void U_inv(AFIELD&, const AFIELD&);
  void Ldag_inv(AFIELD&, const AFIELD&);
  void Udag_inv(AFIELD&, const AFIELD&);

  void mult_gm5_4d(AFIELD&, const AFIELD&);

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Domainwall<AFIELD>(params); }

  static AFopr<AFIELD> *create_object(AFOPR *fopr)
  { return new AFopr_Domainwall<AFIELD>(fopr); }

  static AFopr<AFIELD> *create_object_with_params2(AFOPR *fopr,
                                                   const Parameters& params)
  { return new AFopr_Domainwall<AFIELD>(fopr, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AFOPR::Factory_params::Register(
      "Domainwall", create_object_with_params);

    init &= AFOPR::Factory_fopr::Register(
      "Domainwall", create_object);

    init &= AFOPR::Factory_fopr_params::Register(
      "Domainwall", create_object_with_params2);

    return init;
  }
#endif
};

#endif
