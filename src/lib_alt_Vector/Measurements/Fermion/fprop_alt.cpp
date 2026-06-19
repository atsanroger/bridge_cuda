/*!
        @file    fprop_alt.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib_alt_Vector/Field/afield.h"


// explicit instanciation for AField<double>.

template class Fprop_alt<AField<double,VECTOR> >;

template class Fprop_alt<AField<float,VECTOR> >;

//============================================================END=====
