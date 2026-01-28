/*!
        @file    aforce_F_Smeared.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$
        @version $LastChangedRevision: 2513 $
*/

#ifndef AFORCE_F_SMEARED_INCLUDED
#define AFORCE_F_SMEARED_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Smear/forceSmear.h"
#include "lib/Smear/director_Smear.h"
#include "lib/Fopr/afopr.h"


//! Force calculation for smeared fermion operators.

/*!
    This class determines the force of smeared fermion operator
    using smearing director (MultiSmear instance) and base
    fermion force instance.
                                      [28 Dec 2011 H.Matsufuru]
    Modified: set_mode() is added to incorporate non-hermitian H
                                      [21 Jan 2012 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                      [21 Mar 2015 Y.Namekawa]
*/

template<typename AFIELD>
class AForce_F_Smeared : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;  //!< verbose level

  AForce_F<AFIELD> *m_force;
  ForceSmear *m_force_smear;
  Director_Smear *m_director_smear;
  // Note that this template version always smear the force
  // via Director_Smear.

  using AForce_F<AFIELD>::m_U;
  using AForce_F<AFIELD>::m_Ucp;

 public:
  AForce_F_Smeared(AForce_F<AFIELD> *force,
                   ForceSmear *force_smear,
                   Director_Smear *director_smear)
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

#ifdef USE_FACTORY
 private:
  static AForce_F<AFIELD> *create_object_with_force_forcesmear_director(
                                           AForce_F<AFIELD> *fopr,
                                           ForceSmear* forcesmear,
                                           Director *director)
  { return new AForce_F_Smeared(fopr, forcesmear,
                                (Director_Smear*)director); }

 public:
  static bool register_factory()
  {
    bool init1 = AForce_F<AFIELD>::Factory_force_forcesmear_director::
                  Register("Smeared",
                           create_object_with_force_forcesmear_director);
    return init1;
  }

#endif

};
#endif
