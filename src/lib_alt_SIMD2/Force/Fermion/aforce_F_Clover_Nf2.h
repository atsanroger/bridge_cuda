/*!
        @file    aforce_F_Clover_Nf2.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef AFORCE_F_CLOVER_NF2_INCLUDED
#define AFORCE_F_CLOVER_NF2_INCLUDED

// include files in core core library.
#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Field/field.h"
#include "lib/Parameters/parameters.h"
#include "lib/Fopr/afopr.h"

// include files in alternative code.
#include "lib_alt_SIMD2/Force/Fermion/aforce_F_Wilson_Nf2.h"
#include "lib_alt_SIMD2/Force/Fermion/aforce_F_CloverTerm.h"

//! Force for the standard Clover fermion operator

/*!
    This class calculates the force of the standard Clover
    fermion with two flavors.
    This class is an alternative version of Force_F_Clover_Nf2
    in the core library.
                                     [28 Jan 2019 H.Matusfuru]
 */

template<typename AFIELD>
class AForce_F_Clover_Nf2 : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using AForce_F<AFIELD>::m_vl;
  using AForce_F<AFIELD>::m_U;
  using AForce_F<AFIELD>::m_Ucp;
  static const std::string class_name;

 private:
  // input parameters
  real_t m_kappa;               //!< hopping_parameter
  real_t m_cSW;                 //!< clover_coefficient
  std::vector<int> m_boundary;  //!< boundary_condition
  std::string m_repr;           //!< gamma_matrix_type

  // class objects
  AFopr<AFIELD> *m_fopr;  //!< Clover fermion operator
  AForce_F_Wilson_Nf2<AFIELD> *m_force_w;   //!< Wilson fermion force
  AForce_F_CloverTerm<AFIELD> *m_force_csw; //!< Clover term force

  // working vectors
  AFIELD m_force1;  //!< gauge field
  AFIELD m_force2;  //!< gauge field
  AFIELD m_zeta;    //!< fermion field

 public:
  //! constructor.
  AForce_F_Clover_Nf2(AFopr<AFIELD>* fopr, const Parameters& params)
    : AForce_F<AFIELD>(), m_fopr(fopr)
    { init(params); }

  // destructor.
  ~AForce_F_Clover_Nf2()
    { tidyup(); }

  //! set up parameters with Parameter.
  void set_parameters(const Parameters& params);

  //! set up parameters with values.
  void set_parameters(const double kappa,
                      const double cSW,
                      const std::vector<int> bc);

  //! extracting paramters.
  void extract_parameters(double& kappa,
                          double& cSW,
                          std::vector<int>& bc,
                          const Parameters& params);

  //! setup configuration.
  void set_config(Field *U);

  //! derivative of force wrt. gauge field.
  void force_udiv(AFIELD& force, const AFIELD& eta);

  //! derivative of one component of force wrt. gauge field.
  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:
  //! initial setup.
  void init(const Parameters& params);

  //! final clean-up.
  void tidyup();

  //! implementation of force_udiv1.
  void force_udiv1_impl(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);


#ifdef USE_FACTORY
 private:
  static AForce_F<AFIELD> *create_object_with_fopr_params(
                       AFopr<AFIELD> *fopr, const Parameters& params)
  { return new AForce_F_Clover_Nf2(fopr, params); }

 public:
  static bool register_factory()
  {
    bool init1 = AForce_F<AFIELD>::Factory_fopr_params::Register(
                            "Clover", create_object_with_fopr_params);
    return init1;
  }
#endif

};
#endif
