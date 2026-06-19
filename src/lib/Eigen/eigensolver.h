/*!
        @file    eigensolver.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef EIGENSOLVER_INCLUDED
#define EIGENSOLVER_INCLUDED

//! Eigensolver class for abstract base class of eigen solvers.

/**
   Eigensolver class provides an abstract base class for solvers
   of eigenvalues and eigenvectors.

   This class was changed to template class so as to be available
   with general gield and fermion operator classes.
   The original claas was Eigensolver.
                                       [23 Apr 2018 H.Matsufuru]
 */

#include "Eigen/aeigensolver.h"
#include "Field/field.h"
#include "Fopr/fopr.h"

typedef AEigensolver<Field, Fopr> Eigensolver;

#endif
