/*!
        @file    afopr_Wilson_TwistedMass.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-05 20:48:53 #$

        @version $LastChangedRevision: 2530 $
*/

#ifndef AFOPR_WILSON_TWISTEDMASS_INCLUDED
#define AFOPR_WILSON_TWISTEDMASS_INCLUDED

//#include "fopr_Wilson.h"
#include "Fopr/afopr.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Twisted-mass Wilson fermion operator.

/*!
    This fermion operator defines the twisted-mass Wilson fermion.
    Note that H (= gm5 * D) operator is not hermitian by definition.
                                         [23 Dec 2011 H.Matsufuru]
    YAML is implemented.                 [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                         [21 Mar 2015 Y.Namekawa]
    --> abolished.
    Changed to a template class in ver.2.0.
                                         [19 Dec 2021 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Wilson_TwistedMass : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;

  static const std::string class_name;

 private:
  real_t m_kappa;               //!< hopping parameter
  real_t m_tw_mass;             //!< twisted mass
  std::vector<int> m_boundary;  //!< boundary condition

  std::string m_kernel_type;    //!< kernel type
  std::string m_repr;           //!< gamma matrix representation

  Bridge::VerboseLevel m_vl;    //!< verbose level

  AFopr<AFIELD> *m_fopr_w;      //!< kernel fermion operator

  std::string m_mode;           //!<  mode of multiplication

  // internal data members
  int m_Nc, m_Nd, m_NinF;
  int m_Nvol, m_Ndim;

  const Field_G *m_U;           //!< gauge configuration (pointer)

  AFIELD m_v2;                  //!< working field

  bool m_is_initial_step;       //!< to avoid redundant setup

 public:

  DEPRECATED
  AFopr_Wilson_TwistedMass()
    : m_vl(CommonParameters::Vlevel()) { init("Dirac"); }

  DEPRECATED
  AFopr_Wilson_TwistedMass(const std::string repr) { init(repr); }

  AFopr_Wilson_TwistedMass(const Parameters& params) { init(params); }

  ~AFopr_Wilson_TwistedMass() { tidyup(); }

 private:
  void init(const std::string& repr);

  void init(const Parameters& params);

  void tidyup();

 public:
  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const double tw_mass,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void mult(AFIELD& v, const AFIELD& w);
  void mult_dag(AFIELD& v, const AFIELD& w);

  void set_mode(std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult_gm5(AFIELD& v, const AFIELD& w);

  void mult_gm5p(const int mu, AFIELD& v, const AFIELD& w);

  void D(AFIELD&, const AFIELD&);
  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);
  void Hdag(AFIELD&, const AFIELD&);

  //! returns true if additional field conversion is needed.
  virtual bool needs_convert()
  { return m_fopr_w->needs_convert(); }

  //! converts a Field object into other format if necessary.
  virtual void convert(AFIELD& v, const Field& w)
  { m_fopr_w->convert(v, w); }

  //! reverses to a Field object from other format if necessary.
  virtual void reverse(Field& v, const AFIELD& w)
  { m_fopr_w->reverse(v, w); }

  int field_nin()  { return m_fopr_w->field_nin(); }
  int field_nvol() { return m_fopr_w->field_nvol(); }
  int field_nex()  { return m_fopr_w->field_nex(); }

  //! this returns the number of floating point operations.
  double flop_count();

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object()
  { return new AFopr_Wilson_TwistedMass<AFIELD>(); }

  static AFopr<AFIELD> *create_object_with_repr(const std::string& repr)
  { return new AFopr_Wilson_TwistedMass<AFIELD>(repr); }

  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Wilson_TwistedMass<AFIELD>(params); }

 public:
  static bool register_factory()
  {
    bool                init = true;
    init &= AFopr<AFIELD>::Factory_noarg::Register("Wilson_TwistedMass",
                                                   create_object);
    init &= AFopr<AFIELD>::Factory_string::Register("Wilson_TwistedMass",
                                                    create_object_with_repr);
    init &= AFopr<AFIELD>::Factory_params::Register("Wilson_TwistedMass",
                                                    create_object_with_params);
    return init;
  }
#endif
};
#endif
