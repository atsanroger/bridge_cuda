/*!
      @file    afopr_CloverTetm_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2537 $
 */

#include "lib_alt_Accel/Fopr/afopr_CloverTerm.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Measurements/Gauge/staple_lex.h"
#include "lib/Tools/gammaMatrixSet.h"
#include "lib/Tools/timer.h"
#include "lib/Field/field_G.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/shiftAField_lex.h"
#include "lib_alt_Accel/Measurements/Gauge/astaple_lex.h"

#include "lib_alt/Solver/asolver_CG.h"

typedef float real_t;

// inline function files
#include "lib_alt_Accel/inline/afield_th-inc.h"

#include "lib_alt_Accel/inline/define_params.h"

#include "lib_alt_Accel/Field/afield_Gauge-inc.h"

// function template files
#include "lib_alt_Accel/Field/afield-inc.h"

// library
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Wilson.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Clover.h"

// class template
#include "lib_alt_Accel/Fopr/afopr_CloverTerm-tmpl.h"


template<>
const std::string AFopr_CloverTerm<AField<float,ACCEL> >::class_name
                            = "AFopr_CloverTerm<AField<float,ACCEL> >";

#ifdef USE_FACTORY
namespace {
  //AFopr<AField<float,ACCEL> > *create_object() {
  //return new AFopr_CloverTerm<AField<float,ACCEL> >();
  //}

  //bool init1 = AFopr<AField<float,ACCEL> >::Factory_noarg::Register(
  //                                        "CloverTerm", create_object);

  AFopr<AField<float,ACCEL> > *create_object_with_params(
                                      const Parameters& params) {
    return new AFopr_CloverTerm<AField<float,ACCEL> >(params);
  }

  bool init2 = AFopr<AField<float,ACCEL> >::Factory_params::Register(
                              "CloverTerm", create_object_with_params);
}
#endif


// explicit instanciation.
template class AFopr_CloverTerm<AField<float,ACCEL> >;

//============================================================END=====
