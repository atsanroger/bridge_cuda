/*!
      @file    aprecond.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt/Solver/aprecond.h"

#include "lib_alt_QXS/Field/afield.h"

// explicit instanciation.
template class APrecond<AField<float, QXS> >;
template class APrecond<AField<double, QXS> >;

//============================================================END=====
