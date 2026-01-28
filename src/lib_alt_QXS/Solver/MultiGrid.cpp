/*!

        @file    $Id: MultiGrid.cpp #$

        @brief   base class for MultiGrid

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2499 $

 */

//====================================================================

#include "lib_alt/Solver/MultiGrid.h"
#include "lib_alt_QXS/Field/afield.h"

typedef AField<float, QXS>    AField_f;
typedef AField<double, QXS>   AField_d;

template class MultiGrid<AField_f, AField_f>;
template class MultiGrid<AField_d, AField_d>;

//============================================================END=====
