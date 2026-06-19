/*!
        @file    field_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/


#ifndef FIELD_SF_INCLUDED
#define FIELD_SF_INCLUDED

#include "Parameters/parameters.h"
#include "Field/field_G.h"

//! SU(N) gauge field class in which a few functions are added for the SF.

/*!
  This class defines SU(N) gauge field, which is used such as gauge configuration.
  <ul>
  <li>A derived class from Field_G in order to add a few functions to manipulate the boundary link variables.
  <li>An inheritance was adopted since we need to manipulate the Field contents and upcast into Filed_G object.
  <li>[23 Mar 2012 Y.Taniguchi]
  <li>Field_G_SF and Field_F_SF classes are integrated to Field_SF
     namespace and  converted to functions  [26 Jan 2023 H.Matsufuru]
  </ul>
 */

//----------------------------------------------------------------
// free function version

namespace Field_SF {
  void set_boundary_zero(Field& f);

  void set_boundary_wk(Field_G& u, const Mat_SU_N& wk);

  void set_boundary_wkpr(Field_G& u, const Mat_SU_N& wkpr);

  void set_boundary_zero(Field_G& u);

  void set_boundary_spatial_link_zero(Field_G& u);

  void mult_ct_boundary(Field_G& u, const int t, const double ct);

  void set_boundary_matrix(Mat_SU_N& mat, const std::vector<double>& phi);
}

#endif
