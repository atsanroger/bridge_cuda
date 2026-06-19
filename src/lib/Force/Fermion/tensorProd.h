/*!
        @file    tensorProd.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#ifndef TENSORPROD_INCLUDED
#define TENSORPROD_INCLUDED

#include "Field/field_F.h"

void tensorProd_Field_F(Field_G&, const Field_F&, const Field_F&);
void tensorProd_Field_F(Field_G&, const int mu, const Field_F&, const Field_F&);

#endif
