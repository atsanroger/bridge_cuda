/*!
        @file    fermionFlow_1pt_Function.h

        @brief

        @author  <Yusuke Taniguchi> tanigchi@het.ph.tsukuba.ac.jp
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef FermionFlow_1pt_Function_INCLUDED
#define FermionFlow_1pt_Function_INCLUDED

#include "fprop.h"
//#include "noiseVector_Z2.h"
#include "source_Random.h"
#include "contract_4spinor.h"
#include "corr2pt_4spinor.h"

#include "fermionFlow.h"
#include "Measurements/Gauge/gradientFlow.h"

#include "Tools/gammaMatrixSet.h"
#include "Tools/gammaMatrix.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#include <iostream>
#include <fstream>
#include <sstream>


//! Measures quark one point function with fermion flow.

/*!
  This class measures quark one point function at flow time t.
  The VEV is measured by setting a noise vector at flow time t
  and using adjoint flow for a function "measure_adjoint_flow".
  The noise vector is set at flow time t=0 and uses normal flow
  for a funciton "measure_normal_flow".
  The one point function includes the quark scalar, paseudo-scalar,
  vector and axial vector condensate and those needed to measure
  quark contribution to the energy-momentum tensor.
  At the construction, fermion operator and noise vector generator
  must be specified.
  This is written by Y. Taniguchi and K. Kanaya for ver.1.2.2.
  [05 Jun 2015 Y.Taniguchi]

  For adjoint flow the flowed gauge link is stored at each femion flow step.

t 0---|---|---|---|---|---|---|---|---|
  ^           ^           ^           ^
              ^

 length of "-"=epsilon: m_step_size=0.02
 number of "-" between "|": m_measurement_interval=3
 measurement is perfomed at each "|"
 total number of "|": m_number_of_measurement_times=9
 store gauge conf at each "^"
 number of "|" between "^": m_gauge_store_interval=3
 position of "^^": m_initial_tau=3 for the first measurement
 m_number_of_measurement_times = 9 = last_tau for the last mesurement
 */

class FermionFlow_1pt_Function
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  std::string m_filename_output;

  int m_nc;
  int m_nd;
  int m_ndim;

  std::vector<Fprop *> m_fprop_lex;
  std::vector<Parameters> *m_params_clover;
  Fopr *m_fopr;
  Fprop *m_fprop;
  //  NoiseVector *m_nv;
  Source_Random m_nv;
  //! This object is needed for adjoint flow.
  GradientFlow m_gflow;
  FermionFlow m_fflow;

  //! Number of noise vector for noise estimator.
  int m_Nnoise;
  string m_str_gmset_type;
  double m_hopping_parameter;
  //! Max flow steps = (m_number_of_measurement_times)*(m_measurement_interval)
  int m_number_of_measurement_times;
  //! Steps between adjacent measurements.
  int m_measurement_interval;
  //! Store gauge conf for each t = m_step_size * m_measurement_interval * m_gauge_store_interval
  int m_gauge_store_interval;
  //! start at this tau. tau in physical unit is given by tau=(m_initial_tau)*(m_measurement_interval)*(m_step_size)
  int m_initial_tau;
  //! step size epsilon for Runge-Kutta.
  double m_step_size;

  //! maximum of momentum for Fourier transformation: p_x=[0,max_mom], p_y=[-max_mom,max_mom], p_z=[-max_mom,max_mom]
  int m_max_mom;

  int iex(int emu, int ic, int inoise)
  {
    return emu + (1 + m_ndim) * (ic + m_nc * inoise);
  }

 public:
  FermionFlow_1pt_Function(vector<Fprop *> fprop, Action* action, vector<Parameters> *params_clover)
    : m_vl(CommonParameters::Vlevel()),
    m_nc(CommonParameters::Nc()),
    m_nd(CommonParameters::Nd()),
    m_ndim(CommonParameters::Ndim()),
    m_fprop_lex(fprop),
//    m_nv(nv),
    m_params_clover(params_clover),
    m_gflow(action), m_fflow(action)
  {}

  void set_parameters(const Parameters& params_measurement, const Parameters& params_gflow, const Parameters& params_fflow, const Parameters& params_source_random);
  void set_parameters(const int Nnoise, const int number_of_measurement_times, const int measurement_interval, const int gauge_store_interval, const int initial_tau, const double step_size, const int max_mom);

  void get_parameters(Parameters& params) const;

  //! Measures fermion one point function at flow time t. The VEV is calculated using the standard noise estimator. The estimator is set at flow time t and adjoint flow is used. The same random number is set for each spinor index=0,1,2,3 at one site and color index. The VEV is averaged over four dimensional space-time.
  double measure_adjoint_flow(Field_G& U);

  //! Measures fermion one point function at flow time t. The VEV is calculated using the standard noise estimator. The estimator is set at flow time t=0 and normal fermion flow is used. The same random number is set for each spinor index=0,1,2,3 at one site and color index. The VEV is averaged over four dimensional space-time.
  double measure_normal_order(Field_G& U);

  /*!
    <ul>
    <li>Fermion flow with normal ordering.
    <li>Flow upto t=2.0 for each noise vector.
    <li>Spinor independent source field (same random number is set for each spinor index=0,1,2,3).
    <li>Measure space-time integral of one point functions.
    <li>Measure \vec{x}=0 component of one point functions.
    <li>Measure Fourier transformation of one point functions in space coordinate.
  */
  double measure_disconnected(Field_G& U);

  //! Print out space-time average \f$1/V_{xyzt}\sum_{xyzt}O_i(x,y,z,t)\f$ of quark one point functions for noise vector "inoise" at flow time tt. Inputs are flowed fields xi2, xi and U.
  double print_vev(int inoise, double tt,
                   Field_F_1spinor& xi2, Field_F& xi, Field_G& U);

  //! Print out space average \f$1/V_{xyz}\sum_{xyz}O_i(x,y,z,t)\f$ of quark one point functions as a function of t for noise vector "inoise" at flow time tt. Inputs are flowed fields xi2, xi and U.
  double print_correlator_t(int inoise, double tt,
                            Field_F_1spinor& xi2, Field_F& xi, Field_G& U);

  //! Print out Fourier transformation \f$1/V_{xyz}\sum_{xyz}e^{ip_xx+ip_yy+ip_zz}O_i(x,y,z,t)\f$ of quark one point functions as a function of t for noise vector "inoise" at flow time tt. Inputs are flowed fields xi2, xi and U.
  double print_correlator_t_FT(int inoise, double tt,
                               Field_F_1spinor& xi2, Field_F& xi, Field_G& U);

  //! Print out average \f$1/V_{yzt}\sum_{yzt}O_i(x,y,z,t)\f$ of quark one point functions as a function of t for noise vector "inoise" at flow time tt. Inputs are flowed fields xi2, xi and U.
  double print_correlator_x(int inoise, double tt,
                            Field_F_1spinor& xi2, Field_F& xi, Field_G& U);

  //! Print out Fourier transformation \f$1/V_{yzt}\sum_{yzt}e^{ip_yy+ip_zz+ip_tt}O_i(x,y,z,t)\f$ of quark one point functions as a function of t for noise vector "inoise" at flow time tt. Inputs are flowed fields xi2, xi and U.
  double print_correlator_x_FT(int inoise, double tt,
                               Field_F_1spinor& xi2, Field_F& xi, Field_G& U);

  /*!
   Not implemented yet.
   Print out average \f$1/V_{xzt}\sum_{yzt}O_i(x,y,z,t)\f$ of quark one point functions as a function of t for noise vector "inoise" at flow time tt. Inputs are flowed fields xi2, xi and U.
   */
  double print_correlator_y(int inoise, double tt,
                            Field_F_1spinor& xi2, Field_F& xi, Field_G& U);

  /*!
   Not implemented yet.
   Print out Fourier transformation \f$1/V_{xzt}\sum_{yzt}e^{ip_xx+ip_zz+ip_tt}O_i(x,y,z,t)\f$ of quark one point functions as a function of t for noise vector "inoise" at flow time tt. Inputs are flowed fields xi2, xi and U.
   */
  double print_correlator_y_FT(int inoise, double tt,
                               Field_F_1spinor& xi2, Field_F& xi, Field_G& U);

  /*!
   Not implemented yet.
   Print out average \f$1/V_{xyt}\sum_{yzt}O_i(x,y,z,t)\f$ of quark one point functions as a function of t for noise vector "inoise" at flow time tt. Inputs are flowed fields xi2, xi and U.
   */
  double print_correlator_z(int inoise, double tt,
                            Field_F_1spinor& xi2, Field_F& xi, Field_G& U);

  /*!
   Not implemented yet.
   Print out Fourier transformation \f$1/V_{xyt}\sum_{xyt}e^{ip_xx+ip_yy+ip_tt}O_i(x,y,z,t)\f$ of quark one point functions as a function of t for noise vector "inoise" at flow time tt. Inputs are flowed fields xi2, xi and U.
   */
  double print_correlator_z_FT(int inoise, double tt,
                               Field_F_1spinor& xi2, Field_F& xi, Field_G& U);
};
#endif
