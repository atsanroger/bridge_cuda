/*!
        @file    aforce_F_Smeared_alt.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$
        @version $LastChangedRevision: 2529 $
*/

#ifndef AFORCE_F_SMEARED_ALT_INCLUDED
#define AFORCE_F_SMEARED_ALT_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Smear/aforceSmear.h"
#include "lib_alt/Smear/director_alt_Smear.h"
#include "lib/Fopr/afopr.h"


//! Force calculation for smeared fermion operators.

/*!
    This class determines the force of smeared fermion operator
    using smearing director (MultiSmear instance) and base
    fermion force instance.
                                      [28 Dec 2011 H.Matsufuru]
    AFopr_Smeared_alt receives Director_alt_Smear.
                                      [06 Mar 2023 H.Matsufuru]

*/

template<typename AFIELD>
class AForce_F_Smeared_alt : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  using AForce_F<AFIELD>::m_U;
  static const std::string class_name;

 private:

  AForce_F<AFIELD> *m_force;

  Director_alt_Smear<AFIELD> *m_director_smear;
  // Note that this template version always smear the force
  // via Director_alt_Smear.

  Bridge::VerboseLevel m_vl;  //!< verbose level

  using AForce_F<AFIELD>::m_Ucp;

  AFIELD  m_force1;
  Field_G m_forceG1, m_forceG2;

 public:
  AForce_F_Smeared_alt(AForce_F<AFIELD> *force,
                   Director_alt_Smear<AFIELD> *director_smear)
    : AForce_F<AFIELD>(), m_force(force), m_director_smear(director_smear)
  { init(); }

  void set_parameters(const Parameters&);

  void set_config(Field *U);

  void set_mode(const std::string& mode)
  { m_force->set_mode(mode); }

  void force_udiv(AFIELD& force, const AFIELD& eta);

  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:

  //! initial setup.
  void init();

  void mult_jacobian(Field_G& force);

  void set_config_omp(Field *U);

  void set_config_impl(Field *U);


#ifdef USE_FACTORY
 private:
  static AForce_F<AFIELD> *create_object_with_force_director(
                          AForce_F<AFIELD> *fopr, Director *director)
  { return new AForce_F_Smeared_alt(fopr,
                              (Director_alt_Smear<AFIELD>*)director); }

 public:
  static bool register_factory()
  {
    bool init1 = AForce_F<AFIELD>::Factory_force_director::Register(
                     "Smeared_alt", create_object_with_force_director);
    return init1;
  }

#endif

};
#endif
