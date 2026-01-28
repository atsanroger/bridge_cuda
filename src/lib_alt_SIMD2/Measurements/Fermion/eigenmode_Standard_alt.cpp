/*!
        @file    eigenmode_Standard_alt.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib/ResourceManager/threadManager_OpenMP.h"

#include "lib_alt/Measurements/Fermion/eigenmode_Standard_alt.h"

#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"

#include "lib_alt/Measurements/Fermion/eigenmode_Standard_alt-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD2>.
template<>
const std::string Eigenmode_Standard_alt<AField<double,SIMD2> >::class_name
                         = "Eigenmode_Standard_alt<Afield<double,SIMD2> >";

template class Eigenmode_Standard_alt<AField<double,SIMD2> >;

//====================================================================
// explicit instanciation for AField<float,SIMD2>.
/*
template<>
const std::string Eigenmode_Standard_alt<AField<float,SIMD2> >::class_name
                          = "Eigenmode_Standard_alt<Afield<float,SIMD2> >";

template class Eigenmode_Standard_alt<AField<float,SIMD2> >;
*/

//====================================================================
//============================================================END=====
