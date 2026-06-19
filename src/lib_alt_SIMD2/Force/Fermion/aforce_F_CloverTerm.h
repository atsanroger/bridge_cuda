/*!
        @file    aforce_F_CloverTerm.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef AFORCE_F_CLOVERTERM_INCLUDED
#define AFORCE_F_CLOVERTERM_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Field/field.h"
#include "lib/Parameters/parameters.h"
#include "lib/Tools/gammaMatrix.h"

#include "lib_alt_SIMD2/Field/shiftAField_lex.h"
#include "lib_alt_SIMD2/Measurements/Gauge/astaple_lex.h"


//#include "Fopr/afopr_Clover.h"

//! Force for the standard Clover fermion operator

/*!
    This class calculates the force of the Clover term.
    This class is an alternative version of Force_F_CloverTerm
    in the core library.
                                     [28 Jan 2019 H.Matusfuru]
 */

template<typename AFIELD>
class AForce_F_CloverTerm : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using AForce_F<AFIELD>::m_vl;
  using AForce_F<AFIELD>::m_U;
  using AForce_F<AFIELD>::m_Ucp;
  static const std::string class_name;

  //  AFopr_Clover<AFIELD> *m_fopr_w;
  AFopr<AFIELD> *m_fopr;

 private:
  real_t m_kappa;  //!< hopping parameter
  real_t m_cSW;    //!< clover coefficient
  std::vector<int> m_boundary;  //!< boundary conditions
  std::string m_repr;

  ShiftAField_lex<AFIELD>* m_shiftF;
  ShiftAField_lex<AFIELD>* m_shiftG;
  AStaple_lex<AFIELD>* m_staple;

  int m_Ndim;
  std::vector<GammaMatrix> m_SG;

  AFIELD m_Cud;  //!< working gauge field

  AFIELD m_force1;  //!< working gauge field
  AFIELD m_zeta;    //!< working fermion field
  //  AFIELD m_eta2;    //!< working fermion field
  //  AFIELD m_eta3;    //!< working fermion field

 public:
  AForce_F_CloverTerm(AFopr<AFIELD>* fopr, const Parameters& params)
    : AForce_F<AFIELD>(), m_fopr(fopr)
    { init(params); }

  ~AForce_F_CloverTerm()
    { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa,
                      const double cSW,
                      const std::vector<int> bc);

  void set_config(Field *U);

  void force_udiv(AFIELD& force, const AFIELD& eta);

  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:
  void init(const Parameters& params);

  void tidyup();

  void force_udiv1_impl(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

  void set_component();

  void mult_isigma(AFIELD& v2, const AFIELD& v1, int mu, int nu);

  int index_dir(const int mu, const int nu)
  { return mu + m_Ndim * nu; }

#ifdef USE_FACTORY

  // This class is assumed not to be constructed via Factory.
  /*
 private:
  static AForce_F<AFIELD> *create_object_with_fopr_params(
                       AFopr<AFIELD> *fopr, const Parameters& params)
  { return new AForce_F_CloverTerm(fopr, params); }

 public:
  static bool register_factory()
  {
    bool init1 = AForce_F<AFIELD>::Factory_fopr_params::Register(
                            "CloverTerm", create_object_with_fopr_params);
    return init1;
  }
  */
#endif

};
#endif
