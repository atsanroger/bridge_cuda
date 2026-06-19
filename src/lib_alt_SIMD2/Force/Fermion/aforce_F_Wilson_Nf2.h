/*!
        @file    aforce_F_Wilson_Nf2.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef AFORCE_F_WILSON_NF2_INCLUDED
#define AFORCE_F_WILSON_NF2_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Field/field.h"
#include "lib/Parameters/parameters.h"

#include "lib_alt_SIMD2/Fopr/afopr_Wilson.h"

//! Force for the standard Wilson fermion operator

/*!
    This class calculates the force of the standard Wilson
    fermion with two flavors.
    This class is an alternative version of Force_F_Wilson_Nf2
    in the core library.
                                     [28 Jan 2019 H.Matusfuru]
 */

template<typename AFIELD>
class AForce_F_Wilson_Nf2 : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using AForce_F<AFIELD>::m_vl;
  using AForce_F<AFIELD>::m_U;
  using AForce_F<AFIELD>::m_Ucp;
  static const std::string class_name;

  //  AFopr_Wilson<AFIELD> *m_fopr_w;
  AFopr<AFIELD> *m_fopr_w;

 private:
  double m_kappa;
  std::vector<int> m_boundary;

  AFIELD m_force1;  //!< working gauge field
  AFIELD m_zeta;    //!< working fermion field
  AFIELD m_eta2;    //!< working fermion field
  AFIELD m_eta3;    //!< working fermion field

  std::string m_repr;

 public:
  AForce_F_Wilson_Nf2(AFopr<AFIELD>* fopr, const Parameters& params)
    : AForce_F<AFIELD>(), m_fopr_w(fopr)
    { init(params); }

  ~AForce_F_Wilson_Nf2()
    { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa, const std::vector<int> bc);

  void set_config(Field *U);

  void force_udiv(AFIELD& force, const AFIELD& eta);

  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:
  void init(const Parameters& params);

  void tidyup();

  void force_udiv1_impl(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);


#ifdef USE_FACTORY
 private:
  static AForce_F<AFIELD> *create_object_with_fopr_params(
                       AFopr<AFIELD> *fopr, const Parameters& params)
  { return new AForce_F_Wilson_Nf2(fopr, params); }

 public:
  static bool register_factory()
  {
    bool init1 = AForce_F<AFIELD>::Factory_fopr_params::Register(
                            "Wilson", create_object_with_fopr_params);
    return init1;
  }
#endif

};
#endif
