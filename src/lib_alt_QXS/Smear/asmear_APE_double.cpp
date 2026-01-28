/*!
      @file    asmear_APE_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2513 $
*/

#include "lib_alt_QXS/Smear/asmear_APE.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;


#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

// vector length
#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = ASmear<AField<double,QXS> >::Factory_params::Register(
                                   "APE", create_object_with_params);
}
#endif


template<>
const std::string ASmear_APE<AField<double,QXS> >::class_name
                       = "ASmear_APE<AField<double,QXS> >";


#include "lib_alt_QXS/Smear/asmear_APE-tmpl.h"

// class instanciation.
template class ASmear_APE<AField<double,QXS> >;

//============================================================END=====
