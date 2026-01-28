/*!
      @file    asolver.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Solver/asolver.h"
#include "lib_alt_Accel/Field/afield.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib_alt/Solver/asolver_BiCGStab.h"
#include "lib_alt/Solver/asolver_BiCGStab_Cmplx.h"
#include "lib_alt/Solver/asolver_CG.h"
#include "lib_alt/Solver/asolver_CGNR.h"

template<typename AFIELD>
bool ASolver<AFIELD>::init_factory()
{
  bool result = true;
  result &= ASolver_BiCGStab<AFIELD>::register_factory();
  result &= ASolver_BiCGStab_Cmplx<AFIELD>::register_factory();
  result &= ASolver_CG<AFIELD>::register_factory();
  result &= ASolver_CGNR<AFIELD>::register_factory();
  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


// explicit instanciation.
template class ASolver<AField<float,ACCEL> >;
template class ASolver<AField<double,ACCEL> >;

//============================================================END=====
