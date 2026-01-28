/*!
        @file    aeigensolver.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib/Eigen/aeigensolver.h"

#include "lib/Fopr/afopr.h"
#include "lib_alt_Accel/Field/afield.h"

typedef AField<double,ACCEL> AFIELD;
typedef AFopr<AField<double,ACCEL> > AFOPR;

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib/Eigen/aeigensolver_IRLanczos.h"
#include "lib/Eigen/aeigensolver_IRArnoldi.h"

template<>
bool AEigensolver<AFIELD,AFOPR>::init_factory()
{
  bool result = true;
  result &= AEigensolver_IRLanczos<AFIELD,AFOPR>::register_factory();
  result &= AEigensolver_IRArnoldi<AFIELD,AFOPR>::register_factory();
  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


// explicit instanciation.
template class AEigensolver<AFIELD,AFOPR>;

//============================================================END=====         
