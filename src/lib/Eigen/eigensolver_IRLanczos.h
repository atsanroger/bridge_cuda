/*!
        @file    eigensolver_IRLanczos.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2021-06-15 22:41:26 #$

        @version $LastChangedRevision: 2271 $
*/

#ifndef EIGENSOLVER_IRLANCZOS_INCLUDED
#define EIGENSOLVER_IRLANCZOS_INCLUDED

//! Eigenvalue solver with Implicitly Restarted Lanczos algorithm.

/*!
    This class determines eigenvalues and eigenvectors for a given
    fermion operator.
    Low- or high-lying eigenmodes are determined by chaning
    SortField class object.
                                 [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.         [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                 [21 Mar 2015 Y.Namekawa]
 */

#include "Eigen/aeigensolver_IRLanczos.h"
#include "Eigen/eigensolver.h"
#include "Field/field.h"
#include "Fopr/fopr.h"

typedef AEigensolver_IRLanczos<Field, Fopr> Eigensolver_IRLanczos;

#endif
