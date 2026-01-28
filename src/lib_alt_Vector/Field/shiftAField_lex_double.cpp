/*!
      @file    shiftAField_lex_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/


#include "Field/shiftAField_lex.h"

#include <string>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Vector/Field/aindex_lex.h"
#include "lib_alt_Vector/Field/afield.h"

typedef double real_t;

#include "lib_alt_Vector/inline/define_params.h"

// function implementation
#include "lib_alt_Vector/BridgeVec/bridgeVec_afield.h"
#include "lib_alt_Vector/BridgeVec/bridgeVec_shiftAField_lex.h"

#include "lib_alt_Vector/Field/shiftAField_lex-tmpl.h"

template<>
const std::string ShiftAField_lex<AField<double, VECTOR> >::
         class_name = "ShiftAField_lex<AField<double, VECTOR> >";

// explicit instanciation.
template class ShiftAField_lex<AField<double, VECTOR> >;

//============================================================END=====
