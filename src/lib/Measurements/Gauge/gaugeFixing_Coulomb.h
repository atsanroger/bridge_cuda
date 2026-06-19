/*!
        @file    gaugeFixing_Coulomb.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef GAUGEFIXING_COULOMB_INCLUDED
#define GAUGEFIXING_COULOMB_INCLUDED

#include "gaugeFixing.h"

#include "Field/shiftField_eo.h"
#include "Tools/randomNumbers.h"
#include "Tools/randomNumberManager.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Coulomb gauge fixing.

/*
    This class fixes the gauge of configuration to the Coulomb gauge.
    The implementation assumes that the dimension is 4 and the
    Coulomb gauge fixing is performed within each time slice.
    The algorithm is that developed by the Los Alamos group [see the
    implementation note].
    Overrelaxation is incorporated.
    To escape the Gribov copy, if convergence is not reached on some
    timeslices within the iterations specified by Nreset, random
    gauge transformation is performed to reset the configuration on
    that timeslice.
    This is the reason that random number generator is needed at the
    construction of this class.

    The implementation is not complete:
    - only applies to SU(3) case: because of specific implementation
      of maxTr function (Cabibbo-Marinari maximization).
    - unnecessary arithmetic operations exist for the timeslices
      on which the gauge is already fixed to good precision.
    These should be improved in the version beyond test phase.
                                          [16 Feb 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    Implement YAML.                       [14 Nov 2012 Y.Namekawa]
    Introduce unique_ptr to avoid memory leaks.
                                          [21 Mar 2015 Y.Namekawa]
    Move Staple and RandomNumbers into gaugeFixing.
                                          [30 Mar 2016 Y.Namekawa]
    Add Nc check for USE_GROUP_SU_N.      [31 May 2021 Y.Namekawa]
*/

// strict check
#define CHECK_NC_3                                                 \
  do {                                                             \
    if (CommonParameters::Nc() != 3) {                             \
      vout.crucial(m_vl,                                           \
                   "Error at %s: Nc = 3 is needed, but Nc = %d\n", \
                   class_name.c_str(), CommonParameters::Nc());    \
      exit(EXIT_FAILURE);                                          \
    }                                                              \
  } while (0)

class GaugeFixing_Coulomb : public GaugeFixing
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Niter;             // max iteration number
  int m_Nnaive;            // number of naive iterations
  int m_Nmeas;             // interval of measurements
  int m_Nreset;            // Number of iteration to reset the config.
  double m_Enorm;          // convergence criterion
  double m_wp;             // overrelaxation parameter

  RandomNumbers *m_rand;
  Index_eo m_index;

  ShiftField_eo m_shift;

 public:
  GaugeFixing_Coulomb()
    : m_vl(CommonParameters::Vlevel()),
    m_rand(RandomNumberManager::getInstance())
  {
    CHECK_NC_3;
  }

  GaugeFixing_Coulomb(RandomNumbers *rand)
    : m_vl(CommonParameters::Vlevel()),
    m_rand(rand)
  {
    CHECK_NC_3;
  }

  GaugeFixing_Coulomb(const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_rand(RandomNumberManager::getInstance())
  {
    CHECK_NC_3;

    set_parameters(params);
  }

  GaugeFixing_Coulomb(RandomNumbers *rand, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_rand(rand)
  {
    CHECK_NC_3;

    set_parameters(params);
  }

  ~GaugeFixing_Coulomb() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const int Nnaive,
                      const int Nmeas, const int Nreset,
                      const double Enorm, const double wp);

  void get_parameters(Parameters& params) const;

  void fix(Field_G& Ufix, const Field_G& Uorg);

 private:
  //! one step of gauge fixing with overrelaxation parameter wp.
  void gfix_step(Field_G& Ue, Field_G& Uo, const double wp);

  void set_randomGaugeTrans(const std::valarray<double>& sg, Field_G& Geo);
  void gauge_trans_eo(Field_G& Ue, Field_G& Uo,
                      const Field_G& Geo, const int Ieo);

  void calc_SG(std::valarray<double>& sg, std::valarray<double>& Fval,
               const Field_G& Ue, const Field_G& Uo);
  void calc_DLT(Field_G& Weo,
                const Field_G& Ue, const Field_G& Uo, const int Ieo);
  void calc_W(Field_G& Weo,
              const Field_G& Ue, const Field_G& Uo, const int Ieo);

  void maxTr(Field_G&, Field_G&);
  void maxTr1(Field_G&, Field_G&);
  void maxTr2(Field_G&, Field_G&);
  void maxTr3(Field_G&, Field_G&);

  void sum_global_t(std::valarray<double>& val_global,
                    const std::valarray<double>& val_local);

#ifdef USE_FACTORY
 private:
  static GaugeFixing *create_object()
  {
    return new GaugeFixing_Coulomb();
  }

  static GaugeFixing *create_object_with_params(const Parameters& params)
  {
    return new GaugeFixing_Coulomb(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= GaugeFixing::Factory::Register("Coulomb", create_object);
    init &= GaugeFixing::Factory_params::Register("Coulomb", create_object_with_params);
    return init;
  }
#endif
};
#endif
