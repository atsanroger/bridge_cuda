/*!
      @file    shiftAField_lex_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/


#include "lib_alt_QXS/Field/shiftAField_lex.h"

#include <string>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield.h"

// inline function files
#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/afield_th-inc.h"

// function template files
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt_QXS/Field/shiftAField_lex-tmpl.h"

template<>
const std::string ShiftAField_lex<AField<double, QXS> >::class_name
  = "ShiftAField_lex<AField<double,QXS> >";

// explicit instanciation.
template class ShiftAField_lex<AField<double, QXS> >;

//============================================================END=====
