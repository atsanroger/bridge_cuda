/*!
        @file    gaugeFixing_None.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef GAUGEFIXING_NONE_INCLUDED
#define GAUGEFIXING_NONE_INCLUDED

#include "gaugeFixing.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! None for gauge fixing.

/*
    This class deals with no gauge fixing to manage spectroscopy
    w/o gauge fixing in a single code, proposed by Aoyama-san.
                                        [30 Jun 2016 Y.Namekawa]
*/

class GaugeFixing_None : public GaugeFixing
{
 public:
  static const std::string class_name;

 public:
  GaugeFixing_None() : m_vl(CommonParameters::Vlevel()) {}

  GaugeFixing_None(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

  ~GaugeFixing_None() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const int Nnaive,
                      const int Nmeas, const int Nreset,
                      const double Enorm, const double wp);

  void get_parameters(Parameters& params) const;

  void fix(Field_G& Ufix, const Field_G& Uorg);

 private:
  Bridge::VerboseLevel m_vl;


#ifdef USE_FACTORY
 private:
  static GaugeFixing *create_object()
  {
    return new GaugeFixing_None();
  }

  static GaugeFixing *create_object_with_params(const Parameters& params)
  {
    return new GaugeFixing_None(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= GaugeFixing::Factory::Register("None", create_object);
    init &= GaugeFixing::Factory_params::Register("None", create_object_with_params);
    return init;
  }
#endif
};
#endif
