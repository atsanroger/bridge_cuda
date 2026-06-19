/*!
        @file    eigensolver_IRArnoldi.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Eigen/aeigensolver_IRArnoldi.h"
#include "Eigen/aeigensolver_IRArnoldi-tmpl.h"
#include "Eigen/eigensolver_IRArnoldi.h"

//#include "Field/field.h"
//#include "Fopr/fopr.h"

#define COMPLEX    dcomplex

// explicit instanciation for AField<double>.
template<>
const std::string AEigensolver_IRArnoldi<Field, Fopr>::class_name
  = "Eigensolver_IRArnoldi";

template class AEigensolver_IRArnoldi<Field, Fopr>;


//============================================================END=====
