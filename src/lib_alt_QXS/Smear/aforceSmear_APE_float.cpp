/*!
      @file    aforceSmear_APE_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2513 $
*/

#include "lib_alt_QXS/Smear/aforceSmear_APE.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;


#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

// vector length
#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AProjection<AField<float,QXS> >::Factory_params::Register(
                    "Projection_Stout_SU3", create_object_with_params);
}
#endif


template<>
const std::string AForceSmear_APE<AField<float,QXS> >::class_name
                       = "AforceSmear_APE<AField<float,QXS> >";


#include "lib_alt_QXS/Smear/aforceSmear_APE-tmpl.h"

// class instanciation.
template class AForceSmear_APE<AField<float,QXS> >;

//============================================================END=====
