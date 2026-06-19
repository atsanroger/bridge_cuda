/*!
      @file    astaple_lex_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt_QXS/Measurements/Gauge/astaple_lex.h"

#include <string>

#include "lib/ResourceManager/threadManager.h"


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
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"


// class template
#include "lib_alt_QXS/Measurements/Gauge/astaple_lex-tmpl.h"

template<>
const std::string AStaple_lex<AField<double,QXS> >::class_name
                              = "AStaple_lex<AField<double,QXS> >";

// explicit instanciation.
template class AStaple_lex<AField<double,QXS> >;

//============================================================END=====
