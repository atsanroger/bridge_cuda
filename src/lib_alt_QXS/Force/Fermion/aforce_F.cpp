/*!
        @file    aforce_F.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-11 19:12:24 #$
        @version $LastChangedRevision: 2545 $
*/

#include "lib/Force/Fermion/aforce_F.h"

#include "lib/ResourceManager/threadManager.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib_alt/Force/Fermion/aforce_F_Wilson_Nf2.h"
//#include "lib_alt_QXS/Force/Fermion/aforce_F_Clover_Nf2.h"
//#include "lib_alt/Force/Fermion/aforce_F_Smeared.h"
//#include "lib/Force/Fermion/aforce_F_Rational.h"

template<typename AFIELD>
bool AForce_F<AFIELD>::init_factory()
{
  bool result = true;
  // result &= AForce_F_Wilson_Nf2<AFIELD>::register_factory();
  //  result &= AForce_F_Clover_Nf2<AFIELD>::register_factory();
  //  result &= AForce_F_Smeared<AFIELD>::register_factory();
  //  result &= AForce_F_Rational<AFIELD>::register_factory();
  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */



// include files in alt-code dorectories
#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"

// include files in alt-code dorectories
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"
namespace Alt_Gauge = QXS_Gauge;

typedef AField<double, QXS> AFIELD;

//====================================================================
template<>
void AForce_F<AFIELD>::init()
{
  int Nc   = CommonParameters::Nc();
  int Nin  = 2 * Nc * Nc;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_Ucp.reset(Nin, Nvol, Ndim);
  m_ut.reset( Nin, Nvol, 1);
}

//====================================================================
template<>
void AForce_F<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
// Note that since mult_generator() is called in force_core() and
// force_core1(), its specialization must be placed before the
// latter functions.
template<>
void AForce_F<AFIELD>::mult_generator(AFIELD& force)
{
  //vout.crucial("mult_generator is called.\n");

  int Ndim = force.nex();

#pragma omp barrier

  for(int mu = 0; mu < Ndim; ++mu){
    copy(m_ut, 0, force, mu);
    Alt_Gauge::mult_Gnn(force, mu, m_Ucp, mu, m_ut, 0);
    Alt_Gauge::at_G(force, mu);
  }
  scal(force, -2.0);

  //real_t ff = force.norm2();
  //vout.crucial("norm2 of force (aforce_F) = %e\n", ff);
  
#pragma omp barrier
}

//====================================================================
template<>
void AForce_F<AFIELD>::force_core(AFIELD& force, const AFIELD& eta)
{
  //vout.crucial("force_core is called.\n");
  force_udiv(force, eta);
  mult_generator(force);
}

//====================================================================
template<>
void AForce_F<AFIELD>::force_core1(AFIELD& force,
                                   const AFIELD& zeta,
                                   const AFIELD& eta)
{
  //vout.crucial("force_core1 is called.\n");
  force_udiv1(force, zeta, eta);
  mult_generator(force);
}

//====================================================================

// explicit instanciation.
template class AForce_F<AField<double,QXS> >;
//template class AForce_F<AField<float,QXS> >;

//============================================================END=====
