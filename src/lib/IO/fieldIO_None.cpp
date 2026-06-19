/*!
        @file    fieldIO_None.cpp

        @brief

        @author  $Tatsumi Aoyama (aoym)
                 LastChangedBy: aoym $

        @date    $LastChangedDate: 2014-10-09 17:20:36 #$

        @version $LastChangedRevision: 2314 $
*/

#include "fieldIO_None.h"

#include "bridgeIO.h"
using Bridge::vout;

const std::string FieldIO_None::class_name = "FieldIO_None";

//====================================================================
void FieldIO_None::read_file(Field& v, const std::string& filename)
{
  vout.detailed(m_vl, "read successful: nothing performed.\n");
}


//====================================================================
void FieldIO_None::write_file(Field& v, const std::string& filename)
{
  vout.detailed(m_vl, "write successful: data discarded.\n");
}


//====================================================================
//============================================================END=====
