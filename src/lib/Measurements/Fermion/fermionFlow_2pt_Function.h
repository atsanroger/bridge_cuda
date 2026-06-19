/*!
        @file    fermionFlow_2pt_Function.h

        @brief

        @author  <Yusuke Taniguchi> tanigchi@het.ph.tsukuba.ac.jp
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef FermionFlow_2pt_Function_INCLUDED
#define FermionFlow_2pt_Function_INCLUDED

#include "Field/index_lex.h"
#include "Field/shiftField_lex.h"

#include "fprop.h"
#include "noiseVector_Z2.h"
#include "contract_4spinor.h"
#include "corr2pt_4spinor.h"
#include "Measurements/Fermion/source_Local.h"

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

class FermionFlow_2pt_Function
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
  //! This object is needed for adjoint flow.
  GradientFlow m_gflow;
  FermionFlow m_fflow;

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

  //! Set "backward derivative source" \f$\left(\delta_{x,y+\hat{\mu}}U_\mu^\dagger(y)_{a,ic}-\delta_{x,y-\hat{\mu}}U_\mu(y-\hat{\mu})_{a,ic}\right)/2\f$ at position y, for color index ic and component mu. The field src is intended to be Field_F_1spinor.
  void set_backward_derivative_source(Field& src, const Field_G& U, const std::vector<int>& source_position, const int ic, const int mu);

 public:
  FermionFlow_2pt_Function(vector<Fprop *> fprop, Action* action, vector<Parameters> *params_clover)
    : m_vl(CommonParameters::Vlevel()),
    m_nc(CommonParameters::Nc()),
    m_nd(CommonParameters::Nd()),
    m_ndim(CommonParameters::Ndim()),
    m_fprop_lex(fprop),
    m_params_clover(params_clover),
    m_gflow(action), m_fflow(action)
  {}

  void set_parameters(const Parameters& params_measurement, const Parameters& params_gflow, const Parameters& params_fflow);
  void set_parameters(const int number_of_measurement_times, const int measurement_interval, const int gauge_store_interval, const int initial_tau, const double step_size, const int max_mom);

  void get_parameters(Parameters& params) const;

  /*!
   <ul>
   <li>Measure connected two point functions <O(t,x)O(t,0)>.
   <li>Measure meson connected two point functions in time and space coordinate.
   <li>A point source is set for the propagator.
   <li>Flow upto t for each noise vector.
   <li>Fermion flow with adjoint and normal ordering.
   <li>Measure Fourier transformation of two point functions in space coordinate.
   */
  double measure_meson_correlator(Field_G& U);

  /*!
   <ul>
   <li>Measure connected two point functions <O_{3,5}(t,x)O_{3,5}(t,0)>.
   <li>Measure EM tensor connected two point functions in time and space coordinate.
   <li>A point and backward derivative source is set for the propagator.
   <li>Flow upto t for each noise vector.
   <li>Fermion flow with adjoint and normal ordering.
   <li>Measure Fourier transformation of two point functions in space coordinate.
   */
  double measure_EMT_correlator(Field_G& U);

  //! Calculate meson correlator in t direction. Channels are SS, PP, VV, AA, A_4P, TT.
  double print_meson_correlator_t(double tt,
                                  const vector<vector<Field_F> >& sq);

  //! Calculate meson correlator in t direction including the tree level improvement term. Channels are SS, PP, VV, AA, A_4P, TT.
  double print_meson_correlator_t(double tt,
                                  const vector<vector<Field_F> >& sq,
                                  const vector<Field_F>& cq);

  //! Calculate meson correlator in t direction with Fourier transformation in space. Channels are SS, PP, VV, AA, A_4P, TT.
  double print_meson_correlator_t_FT(double tt,
                                     const vector<vector<Field_F> >& sq);

  //! Calculate meson correlator in t direction with Fourier transformation in space including the tree level improvement term. Channels are SS, PP, VV, AA, A_4P, TT.
  double print_meson_correlator_t_FT(double tt,
                                     const vector<vector<Field_F> >& sq,
                                     const vector<Field_F>& cq);

  //! Calculate meson correlator in x direction. Channels are SS, PP, VV, AA, A_4P, TT.
  double print_meson_correlator_x(double tt,
                                  const vector<vector<Field_F> >& sq);

  //! Calculate meson correlator in x direction including the tree level improvement term. Channels are SS, PP, VV, AA, A_4P, TT.
  double print_meson_correlator_x(double tt,
                                  const vector<vector<Field_F> >& sq,
                                  const vector<Field_F>& cq);

  //! Calculate meson correlator in x direction. Channels are SS, PP, VV, AA, A_4P, TT.
  double print_meson_correlator_x_FT(double tt,
                                     const vector<vector<Field_F> >& sq);

  //! Calculate meson correlator in x direction including the tree level improvement term. Channels are SS, PP, VV, AA, A_4P, TT.
  double print_meson_correlator_x_FT(double tt,
                                     const vector<vector<Field_F> >& sq,
                                     const vector<Field_F>& cq);

  //! Calculate \f$-2{\rm Re tr}\left(\left(\eta_{\sigma}(x)\right)^\dagger\gamma_5\gamma_\mu D_\nu^x\xi(x)\gamma_\rho\gamma_5\right)\f$
  double print_O3O3a_correlator_t(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);
  double print_O3O3a_correlator_x(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);

  //! Calculate \f$-2{\rm Re tr}\left(\left(\eta_{\sigma}(x)\right)^\dagger\gamma_5\gamma_\mu D_\nu^x\xi(x)\gamma_\rho\gamma_5\right)\f$ with Fourier transformation in summed coordinates.
  double print_O3O3a_correlator_t_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);
  double print_O3O3a_correlator_x_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);

  //! Calculate \f$2{\rm Re tr}\left(\left({D}_\nu^x\eta_{\sigma}(x)\right)^\dagger\gamma_5\gamma_\mu\xi(x)\gamma_\rho\gamma_5\right)\f$
  double print_O3O3b_correlator_t(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);
  double print_O3O3b_correlator_x(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);

  //! Calculate \f$2{\rm Re tr}\left(\left({D}_\nu^x\eta_{\sigma}(x)\right)^\dagger\gamma_5\gamma_\mu\xi(x)\gamma_\rho\gamma_5\right)\f$ with Fourier transformation in summed coordinates.
  double print_O3O3b_correlator_t_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);
  double print_O3O3b_correlator_x_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);

  //! Calculate \f$-2{\rm Re tr}\left(\left(\xi(x) \right)^{\dagger}\gamma_{5}\gamma_\mu D_\nu\xi(x)\gamma_{5}\right)\f$
  double print_O3O5_correlator_t(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);
  double print_O3O5_correlator_x(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);

  //! Calculate \f$-2{\rm Re tr}\left(\left(\xi(x) \right)^{\dagger}\gamma_{5}\gamma_\mu D_\nu\xi(x)\gamma_{5}\right)\f$ with Fourier transformation in summed coordinates.
  double print_O3O5_correlator_t_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);
  double print_O3O5_correlator_x_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq);

  //! Calculate \f$2{\rm Re tr}\left(\left(\xi(x)\right)^{\dagger}\gamma_{5}\eta(x)_\nu\gamma_\mu\gamma_{5}\right)\f$
  double print_O5O3_correlator_t(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq);
  double print_O5O3_correlator_x(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq);

  //! Calculate \f$2{\rm Re tr}\left(\left(\xi(x)\right)^{\dagger}\gamma_{5}\eta(x)_\nu\gamma_\mu\gamma_{5}\right)\f$ with Fourier transformation in summed coordinates.
  double print_O5O3_correlator_t_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq);
  double print_O5O3_correlator_x_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq);
};
#endif
