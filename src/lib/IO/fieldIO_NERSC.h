/*!
        @file    fieldIO_NERSC.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2014-04-12 16:31:14 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FIELDIO_NERSC_INCLUDED
#define FIELDIO_NERSC_INCLUDED

#include <string>
using std::string;

#include "fieldIO.h"

#include "Field/index_lex.h"
#include "Field/field.h"

//! FieldIO_NERSC class for file I/O of Field data in binary format.

/*!
   Main ingredients of this class were implemented by T.Yoshie:
     cinput, coutput, big_endian, and byte_swap.
   Incorporated to GaugeConfig_Binary class family by H.Matsufuru.


   The file format treated in this class defined in
     Massimio Di Pierro
     "QCDUTILS"
     arXiv:1202.4813 [hep-lat]

   The binary format is the same as FieldIO_Binary class.
   FieldIUO_NERSC class does not provide error checks against the meta data
   encoded in the header.
                                        [11 Apr 2022 I.Kanamori]
 */

class FieldIO_NERSC : public FieldIO
{
 public:
  static const std::string class_name;

 public:
  FieldIO_NERSC(const IO_Format::Format *format) : FieldIO(format)
  {}

  void read_file(Field& v, const std::string& filename);
  void write_file(Field& v, const std::string& filename);
};
#endif
