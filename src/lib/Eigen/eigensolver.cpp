/*!
        @file    eigensolver.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Eigen/aeigensolver.h"
#include "Field/field.h"
#include "Fopr/fopr.h"

typedef Field          AFIELD;
typedef AFopr<Field>   AFOPR;

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib/Eigen/aeigensolver_IRLanczos.h"
#include "lib/Eigen/aeigensolver_IRArnoldi.h"

template<>
bool AEigensolver<AFIELD, AFOPR>::init_factory()
{
  bool result = true;
  result &= AEigensolver_IRLanczos<AFIELD, AFOPR>::register_factory();
  result &= AEigensolver_IRArnoldi<AFIELD, AFOPR>::register_factory();
  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


// explicit instanciation.
template class AEigensolver<Field, Fopr>;

//typedef AEigensolver<Field> Eigensolver;

//============================================================END=====
