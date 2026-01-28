/*!
        @file    bridge_init_factory_alt_Accel.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
        @version $LastChangedRevision: 2535 $
*/

#include "lib_alt_Accel/bridge_init_factory_alt_Accel.h"

#ifdef USE_FACTORY

// alt-code
#include "lib/Fopr/afopr.h"
#include "lib/Eigen/aeigensolver.h"
//#include "lib_alt_Accel/Force/Fermion/aforce_F.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt/Solver/asolver.h"

// adding to corelib
//#include "lib_alt_Accel/Action/action.h"

/*
#include "Force/Gauge/force_G.h"
#include "Smear/forceSmear.h"
#include "Smear/smear.h"
#include "Smear/projection.h"
#include "Measurements/Gauge/gaugeFixing.h"
#include "Measurements/Gauge/staple.h"
#include "Measurements/Fermion/source.h"
#include "Tools/randomNumbers.h"
#include "Tools/gammaMatrixSet.h"
*/

#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory_alt_Accel()
{
  bool result = true;

  result &= AFopr<AField<double, ACCEL> >::init_factory();
  result &= ASolver<AField<double, ACCEL> >::init_factory();
  // result &= AForce_F<AField<double> >::init_factory();
  result &= AEigensolver<AField<double, ACCEL>,
                  AFopr<AField<double, ACCEL> > >::init_factory();

  result &= AFopr<AField<float, ACCEL> >::init_factory();
  result &= ASolver<AField<float, ACCEL> >::init_factory();

  /*
  result &= Force_G::init_factory();
  result &= ForceSmear::init_factory();
  result &= Solver::init_factory();
  result &= Action::init_factory();
  result &= Projection::init_factory();
  result &= Smear::init_factory();
  result &= GaugeFixing::init_factory();
  result &= Source::init_factory();
  result &= Staple::init_factory();
  result &= RandomNumbers::init_factory();
  result &= GammaMatrixSet::init_factory();
#ifdef USE_FFTWLIB
  result &= FFT::init_factory();
#endif    
  */

  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

void bridge_report_factory_alt_Accel()
{
  vout.general("------------------------------------------------\n");
  vout.general("Alt_Accel: Factory entries\n");
  vout.general("------------------------------------------------\n");

  typedef AFopr<AField<double, ACCEL> > AFOPR_d;
  typedef AFopr<AField<float, ACCEL> > AFOPR_f;

  //AFopr<AField<double, ACCEL> >::Factory_params::print("AFopr(Parameters&)");
  AFOPR_d::Factory_noarg::print("AFopr<double, ACCEL> >(void)");
  AFOPR_d::Factory_params::print("AFopr<double, ACCEL>(Parameters&)");
  AFOPR_d::Factory_fopr::print("AFopr<double, ACCEL>(Fopr*)");
  AFOPR_d::Factory_fopr_director::print("AFopr<double, ACCEL>(Fopr*, Director*)");

  AFOPR_f::Factory_noarg::print("AFopr<float, ACCEL> >(void)");
  AFOPR_f::Factory_params::print("AFopr<float, ACCEL>(Parameters&)");
  AFOPR_f::Factory_fopr::print("AFopr<float, ACCEL>(Fopr*)");
  AFOPR_f::Factory_fopr_director::print("AFopr<float, ACCEL>(Fopr*, Director*)");

  typedef AEigensolver<AField<double, ACCEL>, AFOPR_d> AEIGENSOLVER_d;

   AEIGENSOLVER_d::Factory_fopr::print(
               "AEigensolver<AField<float, ACCEL> >, Afopr>(Fopr*)");

  /*
  Fopr::Factory_noarg::print("Fopr(void)");
  Fopr::Factory_fopr::print("Fopr(Fopr*)");
  Fopr::Factory_fopr_director::print("Fopr(Fopr*, Director*)");
  Fopr::Factory_string::print("Fopr(string)");

  Force_G::Factory::print("Force_G");

  ForceSmear::Factory::print("ForceSmear");

  Solver::Factory::print("Solver");

  Action::Factory::print("Action");

  Projection::Factory::print("Projection");
  Smear::Factory::print("Smear");

  GaugeFixing::Factory::print("GaugeFixing");

  Source::Factory::print("Source");

  Staple::Factory::print("Staple");

  RandomNumbers::Factory_int::print("RandomNumbers(int)");
  RandomNumbers::Factory_file::print("RandomNumbers(string)");

  GammaMatrixSet::Factory::print("GammaMatrixSet");
  */
  vout.general("------------------------------------------------\n");

  return;
}

#endif /* USE_FACTORY */
