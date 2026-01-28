/*!
        @file    eigensolver_IRLanczos.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Eigen/aeigensolver_IRLanczos.h"
#include "Eigen/aeigensolver_IRLanczos-tmpl.h"
#include "Eigen/eigensolver_IRLanczos.h"

//#include "Field/field.h"
//#include "Fopr/fopr.h"

#define COMPLEX    dcomplex

// explicit instanciation for AField<double>.
template<>
const std::string AEigensolver_IRLanczos<Field, Fopr>::class_name
  = "Eigensolver_IRLanczos";

template class AEigensolver_IRLanczos<Field, Fopr>;


//============================================================END=====
