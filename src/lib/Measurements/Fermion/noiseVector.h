/*!
        @file    noiseVector.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/


#ifndef NOISEVECTOR_INCLUDED
#define NOISEVECTOR_INCLUDED

#include "Field/field.h"

#include "IO/bridgeIO.h"

//! Base class for noise vector generator.

/*!
    This is the base class of noise vector generator for
    trace calculations.
    This class only defines the interface.
                                     [2 Sep 2012 H.Matsufuru]
 */

class NoiseVector
{
 protected:
  Bridge::VerboseLevel m_vl;

 public:
  NoiseVector()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~NoiseVector() {}

 private:
  // non-copyable
  NoiseVector(const NoiseVector&);
  NoiseVector& operator=(const NoiseVector&);

 public:
  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  //! setting a noise vector.
  virtual void set(Field& v) = 0;
};
#endif
