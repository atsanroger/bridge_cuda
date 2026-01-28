/*!
        @file    MultiGrid.cpp
        @brief   base class for MultiGrid
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2229 $
 */

//====================================================================

#include "lib_alt/Solver/MultiGrid.h"
#include "lib_alt_Accel/Field/afield.h"

typedef AField<float, ACCEL> AField_f;
typedef AField<double, ACCEL> AField_d;

template class MultiGrid< AField_f, AField_f >;
template class MultiGrid< AField_d, AField_d >;

//============================================================END=====
