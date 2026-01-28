/*!
      @file    ashiftsolver_CG.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2527 $
*/

#include "lib/Solver/ashiftsolver_CG.h"

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib/Solver/ashiftsolver_CG-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string AShiftsolver_CG<AField<double, QXS>,
                                  AFopr<AField<double, QXS> > >::
                class_name = "AShiftsolver_CG<AField<double,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AShiftsolver_CG<AField<double, QXS>,
                               AFopr<AField<double, QXS> > >::
                                                register_factory();
}
#endif

template class AShiftsolver_CG<AField<double, QXS>,
			       AFopr<AField<double, QXS> > >;

//====================================================================
// explicit instanciation for AField<float,QXS>.
template<>
const std::string AShiftsolver_CG<AField<float, QXS>,
                                  AFopr<AField<float, QXS> > >::
                class_name = "AShiftsolver_CG<AField<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AShiftsolver_CG<AField<float, QXS>,
                               AFopr<AField<float, QXS> > >::
                                                register_factory();
}
#endif

template class AShiftsolver_CG<AField<float, QXS>,
			       AFopr<AField<float, QXS> > >;


//============================================================END=====
