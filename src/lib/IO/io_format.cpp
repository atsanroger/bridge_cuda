/*!
        @file    io_format.cpp

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/
#include "io_format.h"
#include "io_format_gauge.h"

/**
   Several Data layouts are predefined in namespace IO_Format.
 */

namespace IO_Format {
  namespace {
    const Trivial_Format trivial_format_ = Trivial_Format();
    // const Gauge::ILDG_Format  ildg_format_    = Gauge::ILDG_Format();
    // const Gauge::JLQCD_Format jlqcd_format_   = Gauge::JLQCD_Format();
  }

  const Format *Trivial = &trivial_format_;
  // const Format *Gauge::ILDG  = &ildg_format_;
  // const Format *Gauge::JLQCD = &jlqcd_format_;
}
