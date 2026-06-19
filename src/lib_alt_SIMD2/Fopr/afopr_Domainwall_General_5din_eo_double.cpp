/*!
      @file    afopr_Domainwall_General_5din_eo_double.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_SIMD2/Fopr/afopr_Domainwall_General_5din_eo.h"

// C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;

// Bridge++ core library header files
#include "lib/ResourceManager/threadManager_OpenMP.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Tools/math_Sign_Zolotarev_Omega.h"

// vector length
#define  VLEN    VLEND
#define  VLEN2   VLEND2

// extended Bridge++ header files
#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/index_eo_alt.h"

// inline functions
#include "lib_alt_SIMD2/inline/afield_th-inc.h"

#include "lib_alt_SIMD2/inline/vsimd_double-inc.h"
#include "lib_alt_SIMD2/inline/vsimd_common_double-inc.h"
#if defined USE_GROUP_SU3
#include "lib_alt_SIMD2/inline/vsimd_Wilson_SU3_double-inc.h"
#endif

#include "lib_alt_SIMD2/Fopr/mult_common_th_simd2-inc.h"
#include "lib_alt_SIMD2/Fopr/mult_Wilson_simd2_parts-inc.h"

// function template files
#include "lib_alt_SIMD2/Field/afield-inc.h"

// template file
#include "lib_alt_SIMD2/Fopr/afopr_Domainwall_General_5din_eo-tmpl.h"


template<>
const std::string AFopr_Domainwall_General_5din_eo<AField<double,SIMD2> >
  ::class_name = "AFopr_Domainwall_General_5din_eo<AField<double,SIMD2> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double,SIMD2> >::Factory_params::Register(
               "Domainwall_General_5din_eo", create_object_with_params1);
}
#endif


// explicit instanciation
template class AFopr_Domainwall_General_5din_eo<AField<double,SIMD2> >;

//============================================================END=====
