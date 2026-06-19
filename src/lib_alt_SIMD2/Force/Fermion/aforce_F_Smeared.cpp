/*!
        @file    aforce_F_Smeared.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Smeared.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager_OpenMP.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// include files in alt-code dorectories
#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"

#include "lib_alt/Force/Fermion/aforce_F_Smeared-tmpl.h"


// explicit instanciation for AField<double,SIMD2>.
template<>
const std::string AForce_F_Smeared<AField<double,SIMD2> >::class_name
                               = "AForce_F_Smeared<Afield<double,SIMD2> >";


template class AForce_F_Smeared<AField<double,SIMD2> >;

//============================================================END=====
