/*!
        @file    quarkNumberSusceptibility_Wilson.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef QUARKNUMBERSUSCEPTIBILITY_WILSON_INCLUDED
#define QUARKNUMBERSUSCEPTIBILITY_WILSON_INCLUDED

#include "fprop.h"
#include "noiseVector_Z2.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Quark number susceptibility for the Wilson-type fermion.

/*!
    This is class measures the traces which is used to determine
    the quark number susceptibility.
    At the construction, fermion operator and noise vector generator
    must be specified.
                                          [02 Sep 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.                  [06 Jun 2013 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                          [21 Mar 2015 Y.Namekawa]
 */

class QuarkNumberSusceptibility_Wilson
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  Fopr *m_fopr;
  Fprop *m_fprop;
  NoiseVector *m_nv;

  int m_Nnoise;

 public:
  QuarkNumberSusceptibility_Wilson(Fopr *fopr, Fprop *fprop, NoiseVector *nv)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr), m_fprop(fprop), m_nv(nv) {}

  QuarkNumberSusceptibility_Wilson(Fopr *fopr, Fprop *fprop, NoiseVector *nv,
                                   const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr(fopr), m_fprop(fprop), m_nv(nv)
  {
    set_parameters(params);
  }

 private:
  // non-copyable
  QuarkNumberSusceptibility_Wilson(const QuarkNumberSusceptibility_Wilson&);
  QuarkNumberSusceptibility_Wilson& operator=(const QuarkNumberSusceptibility_Wilson&);

 public:
  void set_parameters(const Parameters& params);
  void set_parameters(const int Nnoise);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  void get_parameters(Parameters& params) const;

  //! measure tr1 = Tr[D1*Sq], tr2 = Tr[D2*Sq], tr3 = Tr[D1*Sq*D1*Sq].
  double measure();
};
#endif
