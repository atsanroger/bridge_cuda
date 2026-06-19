/*!
        @file    gaugeFixing_None.cpp

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "gaugeFixing_None.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = GaugeFixing_None::register_factory();
}
#endif

const std::string GaugeFixing_None::class_name = "GaugeFixing_None";

//====================================================================
void GaugeFixing_None::set_parameters(const Parameters& params)
{
  //- No parameters are set.
}


//====================================================================
void GaugeFixing_None::get_parameters(Parameters& params) const
{
  //- No other parameters to get

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void GaugeFixing_None::set_parameters(const int Niter, const int Nnaive,
                                      const int Nmeas, const int Nreset,
                                      const double Enorm, const double wp)
{
  //- No parameters are set.
}


//====================================================================
void GaugeFixing_None::fix(Field_G& Ufix, const Field_G& Uorg)
{
  copy(Ufix, Uorg);  // do nothing for gauge fixing, just copy input.
}


//====================================================================
//============================================================END=====
