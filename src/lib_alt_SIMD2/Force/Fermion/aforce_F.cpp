/*!
        @file    aforce_F.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib/Force/Fermion/aforce_F.h"

#include "lib/ResourceManager/threadManager_OpenMP.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib_alt_SIMD2/Force/Fermion/aforce_F_Wilson_Nf2.h"
#include "lib_alt_SIMD2/Force/Fermion/aforce_F_Clover_Nf2.h"
#include "lib_alt/Force/Fermion/aforce_F_Smeared.h"
#include "lib/Force/Fermion/aforce_F_Rational.h"

template<typename AFIELD>
bool AForce_F<AFIELD>::init_factory()
{
  bool result = true;
  result &= AForce_F_Wilson_Nf2<AFIELD>::register_factory();
  result &= AForce_F_Clover_Nf2<AFIELD>::register_factory();
  result &= AForce_F_Smeared<AFIELD>::register_factory();
  result &= AForce_F_Rational<AFIELD>::register_factory();
  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */



// include files in alt-code dorectories
#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"

#define  VLEN    VLEND
#define  VLEN2   VLEND2

#include "lib_alt_SIMD2/inline/afield_th-inc.h"

#include "lib_alt_SIMD2/inline/vsimd_double-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_double-inc.h"

#if defined USE_GROUP_SU3
#include "lib_alt_SIMD2/inline/vsimd_Wilson_SU3_double-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_SU3_double-inc.h"
#endif

#include "lib_alt_SIMD2/Field/afield-inc.h"
#include "lib_alt_SIMD2/Field/afield_Gauge-inc.h"


//====================================================================
template<typename AFIELD>
void AForce_F<AFIELD>::init()
{
  int Nc   = CommonParameters::Nc();
  int Nin  = 2 * Nc * Nc;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_Ucp.reset(Nin, Nvol, Ndim);

}

//====================================================================
template<typename AFIELD>
void AForce_F<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
// Note that since mult_generator() is called in force_core() and
// force_core1(), its specialization must be placed before the
// latter functions.
template<typename AFIELD>
void AForce_F<AFIELD>::mult_generator(AFIELD& force)
{
  int Nc   = CommonParameters::Nc();
  int Nin  = 2 * Nc * Nc;
  int Nvol = force.nvol();
  int Ndim = force.nex();

  AFIELD force2(Nin, Nvol, Ndim);

#pragma omp parallel
 {
  copy(force2, force);
  #pragma omp barrier

  for(int ex = 0; ex < Ndim; ++ex){
    mult_Gnn(force, ex, m_Ucp, ex, force2, ex);
    at_G(force, ex);
  }
  #pragma omp barrier

  scal(force, -2.0);
 }    

}

//====================================================================
template<typename AFIELD>
void AForce_F<AFIELD>::force_core(AFIELD& force, const AFIELD& eta)
{
  force_udiv(force, eta);
  mult_generator(force);
}

//====================================================================
template<typename AFIELD>
void AForce_F<AFIELD>::force_core1(AFIELD& force,
			const AFIELD& zeta, const AFIELD& eta)
{
  force_udiv1(force, zeta, eta);
  mult_generator(force);
}

//====================================================================

// explicit instanciation.
template class AForce_F<AField<double,SIMD2> >;
//template class AForce_F<AField<float> >;


//============================================================END=====
