/*!
        @file    aeigensolver.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-04-01 12:13:09 #$
        @version $LastChangedRevision: 2595 $
*/

#include "lib/Eigen/aeigensolver.h"

#include "lib/Fopr/afopr.h"
#include "lib_alt_QXS/Field/afield.h"

typedef AField<double, QXS>           AFIELD_d;
typedef AFopr<AField<double, QXS> >   AFOPR_d;

typedef AField<float, QXS>            AFIELD_f;
typedef AFopr<AField<float, QXS> >    AFOPR_f;

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib/Eigen/aeigensolver_IRLanczos.h"
#include "lib/Eigen/aeigensolver_IRArnoldi.h"
#include "lib_alt/Eigen/aeigensolver_IRArnoldi_block.h"

template<>
bool AEigensolver<AFIELD_d, AFOPR_d>::init_factory()
{
  typedef AFIELD_d   AFIELD;
  typedef AFOPR_d    AFOPR;
  bool result = true;
  result &= AEigensolver_IRLanczos<AFIELD, AFOPR>::register_factory();
  result &= AEigensolver_IRArnoldi<AFIELD, AFOPR>::register_factory();
  result &= AEigensolver_IRArnoldi_block<AFIELD, AFOPR>::register_factory();
  return result;
}


template<>
bool AEigensolver<AFIELD_f, AFOPR_f>::init_factory()
{
  typedef AFIELD_f   AFIELD;
  typedef AFOPR_f    AFOPR;
  bool result = true;
  result &= AEigensolver_IRLanczos<AFIELD, AFOPR>::register_factory();
  result &= AEigensolver_IRArnoldi<AFIELD, AFOPR>::register_factory();
  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


// explicit instanciation.
template class AEigensolver<AFIELD_d, AFOPR_d>;
template class AEigensolver<AFIELD_f, AFOPR_f>;

//============================================================END=====
