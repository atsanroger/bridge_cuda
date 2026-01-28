/*!
        @file    eigensolver_IRArnoldi.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

        @version $LastChangedRevision: 2422 $
*/

#ifndef EIGENSOLVER_IRARNORDI_INCLUDED
#define EIGENSOLVER_IRARNORDI_INCLUDED

//! Eigenvalue solver with Implicitly Restarted Arnoldi algorithm.

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

#include "Eigen/aeigensolver_IRArnoldi.h"
#include "Eigen/eigensolver.h"
#include "Field/field.h"
#include "Fopr/fopr.h"

typedef AEigensolver_IRArnoldi<Field, Fopr> Eigensolver_IRArnoldi;

#endif
