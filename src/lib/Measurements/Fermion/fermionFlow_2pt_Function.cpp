/*!
        @file    fermionFlow_2pt_Function.cpp

        @brief

        @author  <Yusuke Taniguchi> tanigchi@het.ph.tsukuba.ac.jp
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#include "fermionFlow_2pt_Function.h"

const std::string FermionFlow_2pt_Function::class_name = "FermionFlow_2pt_Function";

//====================================================================
void FermionFlow_2pt_Function::set_parameters(const Parameters& params_measurement, const Parameters& params_gflow, const Parameters& params_fflow)
{
  m_filename_output = params_measurement.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  std::string vlevel;
  if (!params_measurement.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int    Nnoise;
  int    number_of_measurement_times;
  int    measurement_interval;
  int    gauge_store_interval;
  int    initial_tau;
  int    max_mom;
  double step_size;

  int err = 0;
  err += params_measurement.fetch_int("initial_tau", initial_tau);
  err += params_measurement.fetch_int("number_of_measurement_times", number_of_measurement_times);
  err += params_measurement.fetch_int("measurement_interval", measurement_interval);
  err += params_measurement.fetch_int("gauge_store_interval", gauge_store_interval);
  err += params_measurement.fetch_int("max_momentum", max_mom);
  err += params_fflow.fetch_double("step_size", step_size);

  if (err) {
    vout.crucial(m_vl, "%s: fetch error, input parameter not found.\n", class_name.c_str());
    abort();
  }

  assert(number_of_measurement_times >= 1);
  assert(initial_tau > 0);

  set_parameters(number_of_measurement_times, measurement_interval, gauge_store_interval, initial_tau, step_size, max_mom);
  m_gflow.set_parameters(params_gflow);
  m_fflow.set_parameters(params_fflow);
}


//====================================================================
void FermionFlow_2pt_Function::get_parameters(Parameters& params) const
{
  params.set_int("initial_tau", m_initial_tau);
  params.set_int("number_of_measurement_times", m_number_of_measurement_times);
  params.set_int("measurement_interval", m_measurement_interval);
  params.set_int("gauge_store_interval", m_gauge_store_interval);
  params.set_int("max_momentum", m_max_mom);
  params.set_double("step_size", m_step_size);

  Parameters params_gflow;
  m_gflow.get_parameters(params_gflow);
  params.set_Parameters("gradient_flow", params_gflow);

  Parameters params_fflow;
  m_gflow.get_parameters(params_fflow);
  params.set_Parameters("fermion_flow", params_fflow);

  params.set_string("filename_output", m_filename_output);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  return;
}


//====================================================================
void FermionFlow_2pt_Function::set_parameters(const int number_of_measurement_times, const int measurement_interval, const int gauge_store_interval, const int initial_tau, const double step_size, const int max_mom)
{
  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  initial_tau= %d\n", initial_tau);
  vout.general(m_vl, "  number_of_measurement_times = %d\n", number_of_measurement_times);
  vout.general(m_vl, "  measurement_interval = %d\n", measurement_interval);
  vout.general(m_vl, "  gauge_store_interval = %d\n", gauge_store_interval);
  vout.general(m_vl, "  step_size = %lf\n", step_size);
  vout.general(m_vl, "  max_momentum = %d\n", max_mom);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(initial_tau);
  err += ParameterCheck::non_zero(initial_tau);
  err += ParameterCheck::non_negative(number_of_measurement_times);
  err += ParameterCheck::non_zero(number_of_measurement_times);
  err += ParameterCheck::non_negative(measurement_interval);
  err += ParameterCheck::non_negative(gauge_store_interval);
  err += ParameterCheck::non_zero(gauge_store_interval);
  err += ParameterCheck::non_negative(step_size);
  err += ParameterCheck::non_negative(max_mom);

  err += ParameterCheck::is_satisfied(number_of_measurement_times >= gauge_store_interval);
  err += ParameterCheck::is_satisfied(number_of_measurement_times % gauge_store_interval == 0);
  // assert(number_of_measurement_times >= gauge_store_interval);
  // assert(number_of_measurement_times % gauge_store_interval == 0);

  if (err) {
    vout.crucial(m_vl, "%s: parameter range check failed.\n", class_name.c_str());
    abort();
  }

  //- store values
  m_number_of_measurement_times = number_of_measurement_times;
  m_measurement_interval        = measurement_interval;
  m_gauge_store_interval        = gauge_store_interval;
  m_step_size   = step_size;
  m_initial_tau = initial_tau;
  m_max_mom     = max_mom;
}


//====================================================================
double FermionFlow_2pt_Function::measure_meson_correlator(Field_G& U)
{
  int Lx    = CommonParameters::Lx();
  int Ly    = CommonParameters::Ly();
  int Lz    = CommonParameters::Lz();
  int Lt    = CommonParameters::Lt();
  int Nvol  = CommonParameters::Nvol();
  int Nx    = CommonParameters::Nx();
  int Ny    = CommonParameters::Ny();
  int Nz    = CommonParameters::Nz();
  int Nt    = CommonParameters::Nt();
  int ipex  = Communicator::ipe(0);
  int ipey  = Communicator::ipe(1);
  int ipez  = Communicator::ipe(2);
  int ipet  = Communicator::ipe(3);
  int Nsvol = Nvol / Nt;

  assert(Nsvol == Nx * Ny * Nz);

  double t      = 0; //dummy
  double result = 0.0;

  int     number_of_confs = (m_initial_tau + m_number_of_measurement_times - 1) / m_gauge_store_interval;
  Field_G U0(Nvol, m_ndim *(number_of_confs + 1));
  // Store Gauge field at flow time tau=0
  for (int mu = 0; mu < m_ndim; ++mu) {
    U0.setpart_ex(mu, U, mu);
  }
  for (int i = 0; i < number_of_confs; ++i) {
    // Store Gauge field at flow time tau
    for (int j = 0; j < m_gauge_store_interval * m_measurement_interval; j++) {
      m_gflow.evolve(t, U);
    }
    for (int mu = 0; mu < m_ndim; ++mu) {
      U0.setpart_ex(mu + m_ndim * (i + 1), U, mu);
    }
  }
  vout.general(m_vl, "%s::measure_meson_correlator\n", class_name.c_str());

  int    Nconv;
  double diff;

  int Nquark  = m_fprop_lex.size();
  int Nquark2 = m_params_clover->size();
  assert(Nquark == Nquark2);

  int Nin = 2 * m_nc * m_nd;

  Field_F         xi(Nvol, Nquark *m_nc *m_nd);
  Field_F_1spinor xi2(Nvol, m_nc);

  Field v(Nin, Nvol, 1);
  Field w(Nin, Nvol, 1);

  std::vector<std::vector<Field_F> > sq(Nquark);
  std::vector<Field_F>               cq;
  for (int iq = 0; iq < Nquark; ++iq) {
    sq[iq].resize(m_nc * m_nd);
  }
  cq.resize(m_nc * m_nd);

  //meson correlator at flow time 0.
  for (int mu = 0; mu < m_ndim; ++mu) {
    U.setpart_ex(mu, U0, mu);
  }
  for (int id = 0; id < m_nd; id++) {
    for (int ic = 0; ic < m_nc; ++ic) {
      int idx = ic + m_nc * id;
      v.set(0.0);
      if ((ipex == 0) && (ipey == 0) && (ipez == 0) && (ipet == 0)) {
        v.set(2 * idx, 0, 0, 1.0);
      }
      for (int iq = 0; iq < Nquark; ++iq) {
        // w = M^-1 * v
        m_fprop_lex[iq]->invert_D(w, v, Nconv, diff);
        vout.general(m_vl, " iq=%d c=%d s=%d Nconv=%d diff=%.8e\n", iq, ic, id, Nconv, diff);
        m_hopping_parameter = (*m_params_clover)[iq].get_double("hopping_parameter");
        scal(w, 2.0 * m_hopping_parameter);
        for (int isite = 0; isite < Nvol; isite++) {
          for (int icc = 0; icc < m_nc; ++icc) {
            for (int idd = 0; idd < m_nd; ++idd) {
              sq[iq][idx].set_ri(icc, idd, isite, 0, w.cmp(2 * (icc + m_nc * idd), isite, 0), w.cmp(1 + 2 * (icc + m_nc * idd), isite, 0));
            }
          }
        }
      }
    }
  }
  print_meson_correlator_t(0.0, sq);
  print_meson_correlator_x(0.0, sq);

  for (int measure = m_initial_tau; measure < m_initial_tau + m_number_of_measurement_times; ++measure) {
    double tt = m_step_size * measure * m_measurement_interval;

    //Set point source at flow time t.
    xi2.set(0.0);
    for (int ic = 0; ic < m_nc; ++ic) {
      if ((ipex == 0) && (ipey == 0) && (ipez == 0) && (ipet == 0)) {
        xi2.set_ri(ic, 0, ic, 1.0, 0.0);
      }
    }

    m_fflow.set_parameters(measure * m_measurement_interval);
    m_fflow.evolve2(xi2, U0, measure, m_gauge_store_interval);

    for (int mu = 0; mu < m_ndim; ++mu) {
      U.setpart_ex(mu, U0, mu);
    }
    for (int id = 0; id < m_nd; id++) {
      for (int ic = 0; ic < m_nc; ++ic) {
        v.set(0.0);
        for (int isite = 0; isite < xi2.nvol(); isite++) {
          for (int icc = 0; icc < m_nc; icc++) {
            v.set(2 * (icc + m_nc * id), isite, 0, xi2.cmp_r(icc, isite, ic));
            v.set(1 + 2 * (icc + m_nc * id), isite, 0, xi2.cmp_i(icc, isite, ic));
          }
        }
        for (int iq = 0; iq < Nquark; ++iq) {
          // w = M^-1 * v
          m_fprop_lex[iq]->invert_D(w, v, Nconv, diff);
          vout.general(m_vl, " iq=%d c=%d s=%d Nconv=%d diff=%.8e\n", iq, ic, id, Nconv, diff);
          m_hopping_parameter = (*m_params_clover)[iq].get_double("hopping_parameter");
          scal(w, 2.0 * m_hopping_parameter);
          int lex = iq + Nquark * (ic + m_nc * id);
          for (int isite = 0; isite < Nvol; isite++) {
            for (int idd = 0; idd < m_nd; ++idd) {
              for (int icc = 0; icc < m_nc; ++icc) {
                xi.set_ri(icc, idd, isite, lex, w.cmp(2 * (icc + m_nc * idd), isite, 0), w.cmp(1 + 2 * (icc + m_nc * idd), isite, 0));
              }
            }
          }
        }
      }
    }
    m_fflow.evolve_normal_order(xi2, xi, U);

    for (int iq = 0; iq < Nquark; ++iq) {
      for (int i = 0; i < m_nc * m_nd; ++i) {
        sq[iq][i].set(0.0);
      }
    }
    for (int i = 0; i < m_nc * m_nd; ++i) {
      cq[i].set(0.0);
    }
    for (int id = 0; id < m_nd; id++) {
      for (int ic = 0; ic < m_nc; ++ic) {
        int idx = ic + m_nc * id;
        for (int isite = 0; isite < Nvol; isite++) {
          for (int icc = 0; icc < m_nc; ++icc) {
            cq[idx].set_ri(icc, id, isite, 0, xi2.cmp_r(icc, isite, ic), xi2.cmp_i(icc, isite, ic));
            for (int iq = 0; iq < Nquark; ++iq) {
              int lex = iq + Nquark * idx;
              for (int idd = 0; idd < m_nd; ++idd) {
                sq[iq][idx].set_ri(icc, idd, isite, 0, xi.cmp_r(icc, idd, isite, lex), xi.cmp_i(icc, idd, isite, lex));
              }
            }
          }
        }
      }
    }
    result += print_meson_correlator_t(tt, sq);
    result += print_meson_correlator_t(tt, sq, cq);
    result += print_meson_correlator_x(tt, sq);
    result += print_meson_correlator_x(tt, sq, cq);
  } // measure
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::measure_EMT_correlator(Field_G& U)
{
  int Nvol = CommonParameters::Nvol();

  // Set source position=(0,0,0,0)
  std::vector<int> source_position(m_ndim, 0);
  Source_Local     source;

  source.set_parameters(source_position);

  double t; //dummy
  double result = 0.0;

  int     number_of_confs = (m_initial_tau + m_number_of_measurement_times - 1) / m_gauge_store_interval;
  Field_G U0(Nvol, m_ndim *(number_of_confs + 1));
  // Store Gauge field at flow time tau=0
  for (int mu = 0; mu < m_ndim; ++mu) {
    U0.setpart_ex(mu, U, mu);
  }
  for (int i = 0; i < number_of_confs; ++i) {
    // Store Gauge field at flow time tau
    for (int j = 0; j < m_gauge_store_interval * m_measurement_interval; j++) {
      m_gflow.evolve(t, U);
    }
    for (int mu = 0; mu < m_ndim; ++mu) {
      U0.setpart_ex(mu + m_ndim * (i + 1), U, mu);
    }
  }
  vout.general(m_vl, "%s::measure_EMT_correlator\n", class_name.c_str());

  int    Nconv;
  double diff;

  int Nquark  = m_fprop_lex.size();
  int Nquark2 = m_params_clover->size();
  assert(Nquark == Nquark2);

  int Nin = 2 * m_nc * m_nd;

  // Stores one D^-1(point source) and Ndim D^-1(drivative source).
  Field_F xi(Nvol, Nquark *m_nc *m_nd * (m_ndim + 1));
  // Stores flow kernel for one point source and Ndim drivative sources.
  Field_F_1spinor xi2(Nvol, m_nc * (m_ndim + 1));
  // Stores D_mu D^-1(source) temporaly
  Field_F Dxi(Nvol, Nquark *m_nc *m_nd * (m_ndim + 1));
  // Stores D_mu(flow kernel) temporaly
  Field_F_1spinor Dxi2(Nvol, m_nc * (m_ndim + 1));

  // Stores one D^-1(point source) for temporary use
  Field_F sxi(Nvol, Nquark *m_nc *m_nd);
  // Stores flow kernel for one point source for temporary use
  Field_F_1spinor sxi2(Nvol, m_nc);
  // For temporary working space
  Field_F_1spinor xscr(Nvol, 1);

  Field v(Nin, Nvol, 1);
  Field w(Nin, Nvol, 1);

  // Store flowed propagators for the function meson_corr()
  // <chi bchi>=sq[0][iq][idx], <chi bchi>D_mu=sq[mu>0][iq][idx]
  std::vector<std::vector<std::vector<Field_F> > > sq(m_ndim + 1);
  // D_mu<chi bchi>=dsq[mu][0][iq][idx], D_mu<chi bchi>D_nu=dsq[mu][nu>0][iq][idx]
  std::vector<std::vector<std::vector<std::vector<Field_F> > > > dsq(m_ndim);
  // <chi bchi>=cq[0][idx], <chi bchi>D_mu=cq[mu>0][idx]
  std::vector<std::vector<Field_F> > cq(m_ndim + 1);
  // D_mu<chi bchi>=dcq[mu][0][idx], D_mu<chi bchi>D_nu=dcq[mu][nu>0][idx]
  std::vector<std::vector<std::vector<Field_F> > > dcq(m_ndim);
  for (int mu = 0; mu < m_ndim + 1; mu++) {
    sq[mu].resize(Nquark);
    for (int iq = 0; iq < Nquark; ++iq) {
      sq[mu][iq].resize(m_nc * m_nd);
    }
    cq[mu].resize(m_nc * m_nd);
  }
  for (int mu = 0; mu < m_ndim; mu++) {
    dsq[mu].resize(m_ndim + 1);
    dcq[mu].resize(m_ndim + 1);
    for (int nu = 0; nu < m_ndim + 1; nu++) {
      dsq[mu][nu].resize(Nquark);
      for (int iq = 0; iq < Nquark; ++iq) {
        dsq[mu][nu][iq].resize(m_nc * m_nd);
      }
      dcq[mu][nu].resize(m_nc * m_nd);
    }
  }

  for (int measure = m_initial_tau; measure < m_initial_tau + m_number_of_measurement_times; ++measure) {
    double tt = m_step_size * measure * m_measurement_interval;
    //    xi2.set(0.0);
    // Point source
    for (int ic = 0; ic < m_nc; ++ic) {
      source.set(xscr, ic);
      xi2.setpart_ex(ic, xscr, 0);
    }
    // Get link variable U at flow time (measure)
    for (int mu = 0; mu < m_ndim; ++mu) {
      U.setpart_ex(mu, U0, mu + m_ndim * (measure / m_gauge_store_interval));
    }
    if (measure % m_gauge_store_interval != 0) {
      for (int i = 0; i < (measure % m_gauge_store_interval) * m_measurement_interval; i++) {
        m_gflow.evolve(t, U);
      }
    }
    // Backward derivative source
    for (int mu = 0; mu < m_ndim; ++mu) {
      for (int ic = 0; ic < m_nc; ++ic) {
        set_backward_derivative_source(xscr, U, source_position, ic, mu);
        xi2.setpart_ex(ic + m_nc * (mu + 1), xscr, 0);
      }
    }

    m_fflow.set_parameters(measure * m_measurement_interval);
    m_fflow.evolve2(xi2, U0, measure, m_gauge_store_interval);

    for (int mu = 0; mu < m_ndim + 1; ++mu) {
      for (int nu = 0; nu < m_ndim; ++nu) {
        U.setpart_ex(nu, U0, nu);
      }
      for (int id = 0; id < m_nd; id++) {
        for (int ic = 0; ic < m_nc; ++ic) {
          v.set(0.0);
          for (int isite = 0; isite < xi2.nvol(); isite++) {
            for (int icc = 0; icc < m_nc; icc++) {
              v.set(2 * (icc + m_nc * id), isite, 0, xi2.cmp_r(icc, isite, ic + m_nc * mu));
              v.set(1 + 2 * (icc + m_nc * id), isite, 0, xi2.cmp_i(icc, isite, ic + m_nc * mu));
            }
          }
          for (int iq = 0; iq < Nquark; ++iq) {
            // w = M^-1 * v
            m_fprop_lex[iq]->invert_D(w, v, Nconv, diff);
            vout.general(m_vl, " iq=%d mu=%d c=%d s=%d Nconv=%d diff=%.8e\n", iq, mu, ic, id, Nconv, diff);
            m_hopping_parameter = (*m_params_clover)[iq].get_double("hopping_parameter");
            scal(w, 2.0 * m_hopping_parameter);
            //int lex=iq+Nquark*(ic+m_nc*(id+m_nd*mu));
            int lex = iq + Nquark * (ic + m_nc * id);

            /*
            for(int isite=0; isite<Nvol; isite++){
              for(int idd=0; idd<m_nd; ++idd){
                for(int icc=0; icc<m_nc; ++icc){
                  xi.set_ri(icc, idd, isite, lex, w.cmp(2*(icc+m_nc*idd), isite, 0), w.cmp(1+2*(icc+m_nc*idd), isite, 0));
                }
              }
            }
            */
            sxi.setpart_ex(lex, w, 0);
          }
        }
      }
      for (int ic = 0; ic < m_nc; ++ic) {
        sxi2.setpart_ex(ic, xi2, ic + m_nc * mu);
      }
      m_fflow.evolve_normal_order(sxi2, sxi, U);
      for (int ic = 0; ic < m_nc; ++ic) {
        xi2.setpart_ex(ic + m_nc * mu, sxi2, ic);
        for (int id = 0; id < m_nd; id++) {
          for (int iq = 0; iq < Nquark; ++iq) {
            xi.setpart_ex(iq + Nquark * (ic + m_nc * (id + m_nd * mu)), sxi, iq + Nquark * (ic + m_nc * id));
          }
        }
      }
    }
    //m_fflow.evolve_normal_order(xi2, xi, U);

    // Get link variable U at flow time (measure)
    for (int mu = 0; mu < m_ndim; ++mu) {
      U.setpart_ex(mu, U0, mu + m_ndim * (measure / m_gauge_store_interval));
    }
    if (measure % m_gauge_store_interval != 0) {
      for (int i = 0; i < (measure % m_gauge_store_interval) * m_measurement_interval; i++) {
        m_gflow.evolve(t, U);
      }
    }

    // Initialize cq and dcq since all the spin components are not filled.
    for (int mu = 0; mu < m_ndim + 1; mu++) {
      for (int i = 0; i < m_nc * m_nd; ++i) {
        cq[mu][i].set(0.0);
      }
    }
    for (int mu = 0; mu < m_ndim; mu++) {
      for (int nu = 0; nu < m_ndim + 1; nu++) {
        for (int i = 0; i < m_nc * m_nd; ++i) {
          dcq[mu][nu][i].set(0.0);
        }
      }
    }
    // Set propagators into sq and cq
    for (int mu = 0; mu < m_ndim + 1; mu++) {
      for (int id = 0; id < m_nd; id++) {
        for (int ic = 0; ic < m_nc; ++ic) {
          int idx  = ic + m_nc * id;
          int lexc = ic + m_nc * mu;
          for (int isite = 0; isite < Nvol; isite++) {
            for (int icc = 0; icc < m_nc; ++icc) {
              cq[mu][idx].set_ri(icc, id, isite, 0, xi2.cmp_r(icc, isite, lexc), xi2.cmp_i(icc, isite, lexc));
              for (int iq = 0; iq < Nquark; ++iq) {
                int lexs = iq + Nquark * (ic + m_nc * (id + m_nd * mu));
                for (int idd = 0; idd < m_nd; ++idd) {
                  sq[mu][iq][idx].set_ri(icc, idd, isite, 0, xi.cmp_r(icc, idd, isite, lexs), xi.cmp_i(icc, idd, isite, lexs));
                }
              }
            }
          }
        }
      }
    }
    result += print_meson_correlator_t(tt, sq[0]);
    result += print_meson_correlator_t(tt, sq[0], cq[0]);
    result += print_meson_correlator_x(tt, sq[0]);
    result += print_meson_correlator_x(tt, sq[0], cq[0]);
    result += print_meson_correlator_t_FT(tt, sq[0]);
    result += print_meson_correlator_t_FT(tt, sq[0], cq[0]);
    result += print_meson_correlator_x_FT(tt, sq[0]);
    result += print_meson_correlator_x_FT(tt, sq[0], cq[0]);

    // Set propagators into dsq and dcq
    for (int mu = 0; mu < m_ndim; mu++) {
      m_fflow.del_symmetric(mu, Dxi, xi, U);
      m_fflow.del_symmetric(mu, Dxi2, xi2, U);
      for (int nu = 0; nu < m_ndim + 1; nu++) {
        for (int id = 0; id < m_nd; id++) {
          for (int ic = 0; ic < m_nc; ++ic) {
            int idx  = ic + m_nc * id;
            int lexc = ic + m_nc * nu;
            for (int isite = 0; isite < Nvol; isite++) {
              for (int icc = 0; icc < m_nc; ++icc) {
                dcq[mu][nu][idx].set_ri(icc, id, isite, 0, Dxi2.cmp_r(icc, isite, lexc), Dxi2.cmp_i(icc, isite, lexc));
                for (int iq = 0; iq < Nquark; ++iq) {
                  int lexs = iq + Nquark * (ic + m_nc * (id + m_nd * nu));
                  for (int idd = 0; idd < m_nd; ++idd) {
                    dsq[mu][nu][iq][idx].set_ri(icc, idd, isite, 0, Dxi.cmp_r(icc, idd, isite, lexs), Dxi.cmp_i(icc, idd, isite, lexs));
                  }
                }
              }
            }
          }
        }
      }
    }
    result += print_O3O3a_correlator_t(tt, sq, cq, dsq, dcq);
    result += print_O3O3b_correlator_t(tt, sq, cq, dsq, dcq);
    result += print_O3O5_correlator_t(tt, sq, cq, dsq, dcq);
    result += print_O5O3_correlator_t(tt, sq, cq);
    result += print_O3O3a_correlator_x(tt, sq, cq, dsq, dcq);
    result += print_O3O3b_correlator_x(tt, sq, cq, dsq, dcq);
    result += print_O3O5_correlator_x(tt, sq, cq, dsq, dcq);
    result += print_O5O3_correlator_x(tt, sq, cq);
    result += print_O3O3a_correlator_t_FT(tt, sq, cq, dsq, dcq);
    result += print_O3O3b_correlator_t_FT(tt, sq, cq, dsq, dcq);
    result += print_O3O5_correlator_t_FT(tt, sq, cq, dsq, dcq);
    result += print_O5O3_correlator_t_FT(tt, sq, cq);
    result += print_O3O3a_correlator_x_FT(tt, sq, cq, dsq, dcq);
    result += print_O3O3b_correlator_x_FT(tt, sq, cq, dsq, dcq);
    result += print_O3O5_correlator_x_FT(tt, sq, cq, dsq, dcq);
    result += print_O5O3_correlator_x_FT(tt, sq, cq);
  } // measure
  return result;
}


//====================================================================
void FermionFlow_2pt_Function::set_backward_derivative_source(Field& src, const Field_G& U, const std::vector<int>& source_position, const int ic, const int mu)
{
  src.set(0.0);
  int Nvol = CommonParameters::Nvol();

  Index_lex        l_index;  // lexical only.
  ShiftField_lex   l_shift;
  Field_G          U1(Nvol, m_ndim);
  std::vector<int> l_source_position(m_ndim);

  std::vector<int> Lsize(m_ndim);
  Lsize[0] = CommonParameters::Lx();
  Lsize[1] = CommonParameters::Ly();
  Lsize[2] = CommonParameters::Lz();
  Lsize[3] = CommonParameters::Lt();
  std::vector<int> Nsize(m_ndim);
  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();
  std::vector<int> ipe(m_ndim);
  ipe[0] = Communicator::ipe(0);
  ipe[1] = Communicator::ipe(1);
  ipe[2] = Communicator::ipe(2);
  ipe[3] = Communicator::ipe(3);


  // Set U_mu^dagger(x)_{icc,ic} at point x+mu
  l_source_position      = source_position;
  l_source_position[mu] += 1;
  l_shift.forward(U1, U, mu);
  for (int nu = 0; nu < m_ndim; ++nu) {
    while (l_source_position[nu] < 0)
    {
      l_source_position[nu] += Lsize[nu];
    }
    l_source_position[nu] %= Lsize[nu];
  }
  bool in_node = true;
  for (int i = 0; i < m_ndim; ++i) {
    int inode = l_source_position[i] / Nsize[i];
    if (inode != ipe[i]) {
      in_node = false;
      break;
    }
  }
  if (in_node) {
    for (int i = 0; i < m_ndim; ++i) {
      l_source_position[i] %= Nsize[i];
    }
    const int isite = l_index.site(l_source_position[0], l_source_position[1], l_source_position[2], l_source_position[3]);
    // Set 1/2*U_mu(isite-mu)^dagger_{icc,ic}=1/2*U_mu(isite-mu)^*_{ic,icc}
    for (int icc = 0; icc < m_nc; ++icc) {
      int    cc   = icc + m_nc * ic;
      double real = U1.cmp_r(cc, isite, mu) * 0.5;
      double imag = U1.cmp_i(cc, isite, mu) * (-0.5);
      src.set(2 * icc, isite, 0, real);
      src.set(1 + 2 * icc, isite, 0, imag);
    }
  }

  // Set U_mu(x-mu)_{icc,ic} at point x-mu
  l_source_position      = source_position;
  l_source_position[mu] -= 1;
  for (int nu = 0; nu < m_ndim; ++nu) {
    while (l_source_position[nu] < 0)
    {
      l_source_position[nu] += Lsize[nu];
    }
    l_source_position[nu] %= Lsize[nu];
  }
  in_node = true;
  for (int i = 0; i < m_ndim; ++i) {
    int inode = l_source_position[i] / Nsize[i];
    if (inode != ipe[i]) {
      in_node = false;
      break;
    }
  }
  if (in_node) {
    for (int i = 0; i < m_ndim; ++i) {
      l_source_position[i] %= Nsize[i];
    }
    const int isite = l_index.site(l_source_position[0], l_source_position[1], l_source_position[2], l_source_position[3]);
    // Set -1/2*U_mu(isite)_{icc,ic}
    for (int icc = 0; icc < m_nc; ++icc) {
      int    cc   = ic + m_nc * icc;
      double real = U1.cmp_r(cc, isite, mu) * (-0.5);
      double imag = U1.cmp_i(cc, isite, mu) * (-0.5);
      src.set(2 * icc, isite, 0, real);
      src.set(1 + 2 * icc, isite, 0, imag);
    }
  }
}


//====================================================================
double FermionFlow_2pt_Function::print_meson_correlator_t(double tt,
                                                          const vector<vector<Field_F> >& sq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = sq.size();

  assert(Nquark == 2);
  GammaMatrix *gm    = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *g5mu  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm5   = new GammaMatrix[Nquark];
  GammaMatrix *gm0   = new GammaMatrix[Nquark];
  GammaMatrix *sigma = new GammaMatrix[m_ndim * (m_ndim - 1) / 2 * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i] = gmset->get_GM(gmset->GAMMA5);
    gm0[i] = gmset->get_GM(gmset->UNITY);
    sigma[0 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA12);
    sigma[1 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA23);
    sigma[2 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA31);
    sigma[3 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA41);
    sigma[4 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA42);
    sigma[5 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA43);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator(mcorr, gm0[iq], gm0[iq], sq[iq], sq[iq]);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "SS_tcorr_ss_iq%d%d %lf %4d %0.16e %0.16e\n",
                   iq, iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  corr.meson_correlator(mcorr, gm0[1], gm0[0], sq[1], sq[0]);
  for (int t = 0; t < Lt; ++t) {
    vout.general(m_vl, "SS_tcorr_ss_iq%d%d %lf %4d %0.16e %0.16e\n",
                 0, 1, tt, t, real(mcorr[t]), imag(mcorr[t]));
  }
  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator(mcorr, gm5[iq], gm5[iq], sq[iq], sq[iq]);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "SS_tcorr_pp_iq%d%d %lf %4d %0.16e %0.16e\n",
                   iq, iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  corr.meson_correlator(mcorr, gm5[1], gm5[0], sq[1], sq[0]);
  for (int t = 0; t < Lt; ++t) {
    vout.general(m_vl, "SS_tcorr_pp_iq%d%d %lf %4d %0.16e %0.16e\n",
                 0, 1, tt, t, real(mcorr[t]), imag(mcorr[t]));
  }
  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], sq[iq], sq[iq]);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "SS_tcorr_vv_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                     iq, iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
    corr.meson_correlator(mcorr, gm[mu + m_ndim], gm[mu], sq[1], sq[0]);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "SS_tcorr_vv_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                   0, 1, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator(mcorr, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], sq[iq], sq[iq]);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "SS_tcorr_aa_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                     iq, iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
    corr.meson_correlator(mcorr, g5mu[mu + m_ndim], g5mu[mu], sq[1], sq[0]);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "SS_tcorr_aa_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                   0, 1, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator(mcorr, g5mu[3 + m_ndim * iq], gm5[iq], sq[iq], sq[iq]);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "SS_tcorr_a3p_iq%d%d %lf %4d %0.16e %0.16e\n",
                   iq, iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  corr.meson_correlator(mcorr, g5mu[3 + m_ndim], gm5[0], sq[1], sq[0]);
  for (int t = 0; t < Lt; ++t) {
    vout.general(m_vl, "SS_tcorr_a3p_iq%d%d %lf %4d %0.16e %0.16e\n",
                 0, 1, tt, t, real(mcorr[t]), imag(mcorr[t]));
  }
  for (int mu = 0; mu < m_ndim * (m_ndim - 1) / 2; ++mu) {
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator(mcorr, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sq[iq], sq[iq]);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "SS_tcorr_tt_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                     iq, iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
    corr.meson_correlator(mcorr, sigma[mu + m_ndim * (m_ndim - 1) / 2], sigma[mu], sq[1], sq[0]);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "SS_tcorr_tt_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                   0, 1, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  delete [] sigma;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_meson_correlator_t(double tt,
                                                          const vector<vector<Field_F> >& sq,
                                                          const vector<Field_F>& cq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = sq.size();

  assert(Nquark > 1);
  GammaMatrix *gm    = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *g5mu  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm5   = new GammaMatrix[Nquark];
  GammaMatrix *gm0   = new GammaMatrix[Nquark];
  GammaMatrix *sigma = new GammaMatrix[m_ndim * (m_ndim - 1) / 2 * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i] = gmset->get_GM(gmset->GAMMA5);
    gm0[i] = gmset->get_GM(gmset->UNITY);
    sigma[0 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA12);
    sigma[1 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA23);
    sigma[2 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA31);
    sigma[3 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA41);
    sigma[4 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA42);
    sigma[5 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA43);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  corr.meson_correlator(mcorr, gm0[0], gm0[0], cq, cq);
  for (int t = 0; t < Lt; ++t) {
    vout.general(m_vl, "CC_tcorr %lf %4d %0.16e %0.16e\n",
                 tt, t, real(mcorr[t]), imag(mcorr[t]));
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator(mcorr, gm0[iq], gm0[iq], sq[iq], cq);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "CS_tcorr_ss_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
    corr.meson_correlator(mcorr, gm0[iq], gm0[iq], cq, sq[iq]);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "SC_tcorr_ss_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator(mcorr, gm5[iq], gm5[iq], sq[iq], cq);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "CS_tcorr_pp_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
    corr.meson_correlator(mcorr, gm5[iq], gm5[iq], cq, sq[iq]);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "SC_tcorr_pp_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], sq[iq], cq);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "CS_tcorr_vv_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
      corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], cq, sq[iq]);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "SC_tcorr_vv_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator(mcorr, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], sq[iq], cq);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "CS_tcorr_aa_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
      corr.meson_correlator(mcorr, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], cq, sq[iq]);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "SC_tcorr_aa_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator(mcorr, g5mu[3 + m_ndim * iq], gm5[iq], sq[iq], cq);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "CS_tcorr_a3p_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
    corr.meson_correlator(mcorr, g5mu[3 + m_ndim * iq], gm5[iq], cq, sq[iq]);
    for (int t = 0; t < Lt; ++t) {
      vout.general(m_vl, "SC_tcorr_a3p_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }

  for (int mu = 0; mu < m_ndim * (m_ndim - 1) / 2; ++mu) {
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator(mcorr, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sq[iq], cq);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "CS_tcorr_tt_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
      corr.meson_correlator(mcorr, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], cq, sq[iq]);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "SC_tcorr_tt_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  delete [] sigma;
  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  delete gmset;

  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_meson_correlator_x(double tt,
                                                          const vector<vector<Field_F> >& sq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = sq.size();

  assert(Nquark == 2);
  GammaMatrix *gm    = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *g5mu  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm5   = new GammaMatrix[Nquark];
  GammaMatrix *gm0   = new GammaMatrix[Nquark];
  GammaMatrix *sigma = new GammaMatrix[m_ndim * (m_ndim - 1) / 2 * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i] = gmset->get_GM(gmset->GAMMA5);
    gm0[i] = gmset->get_GM(gmset->UNITY);
    sigma[0 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA12);
    sigma[1 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA23);
    sigma[2 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA31);
    sigma[3 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA41);
    sigma[4 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA42);
    sigma[5 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA43);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator_x(mcorr, gm0[iq], gm0[iq], sq[iq], sq[iq]);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "SS_xcorr_ss_iq%d%d %lf %4d %0.16e %0.16e\n",
                   iq, iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  corr.meson_correlator_x(mcorr, gm0[1], gm0[0], sq[1], sq[0]);
  for (int t = 0; t < Lx; ++t) {
    vout.general(m_vl, "SS_xcorr_ss_iq%d%d %lf %4d %0.16e %0.16e\n",
                 0, 1, tt, t, real(mcorr[t]), imag(mcorr[t]));
  }
  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator_x(mcorr, gm5[iq], gm5[iq], sq[iq], sq[iq]);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "SS_xcorr_pp_iq%d%d %lf %4d %0.16e %0.16e\n",
                   iq, iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  corr.meson_correlator_x(mcorr, gm5[1], gm5[0], sq[1], sq[0]);
  for (int t = 0; t < Lx; ++t) {
    vout.general(m_vl, "SS_xcorr_pp_iq%d%d %lf %4d %0.16e %0.16e\n",
                 0, 1, tt, t, real(mcorr[t]), imag(mcorr[t]));
  }
  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], sq[iq], sq[iq]);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "SS_xcorr_vv_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                     iq, iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
    corr.meson_correlator_x(mcorr, gm[mu + m_ndim], gm[mu], sq[1], sq[0]);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "SS_xcorr_vv_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                   0, 1, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator_x(mcorr, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], sq[iq], sq[iq]);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "SS_xcorr_aa_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                     iq, iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
    corr.meson_correlator_x(mcorr, g5mu[mu + m_ndim], g5mu[mu], sq[1], sq[0]);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "SS_xcorr_aa_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                   0, 1, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator_x(mcorr, g5mu[m_ndim * iq], gm5[iq], sq[iq], sq[iq]);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "SS_xcorr_a0p_iq%d%d %lf %4d %0.16e %0.16e\n",
                   iq, iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  corr.meson_correlator_x(mcorr, g5mu[m_ndim], gm5[0], sq[1], sq[0]);
  for (int t = 0; t < Lx; ++t) {
    vout.general(m_vl, "SS_xcorr_a0p_iq%d%d %lf %4d %0.16e %0.16e\n",
                 0, 1, tt, t, real(mcorr[t]), imag(mcorr[t]));
  }
  for (int mu = 0; mu < m_ndim * (m_ndim - 1) / 2; ++mu) {
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator_x(mcorr, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sq[iq], sq[iq]);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "SS_xcorr_tt_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                     iq, iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
    corr.meson_correlator_x(mcorr, sigma[mu + m_ndim * (m_ndim - 1) / 2], sigma[mu], sq[1], sq[0]);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "SS_xcorr_tt_iq%d%d %d %lf %4d %0.16e %0.16e\n",
                   0, 1, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  delete [] sigma;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_meson_correlator_x(double tt,
                                                          const vector<vector<Field_F> >& sq,
                                                          const vector<Field_F>& cq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = sq.size();

  assert(Nquark > 1);
  GammaMatrix *gm    = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *g5mu  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm5   = new GammaMatrix[Nquark];
  GammaMatrix *gm0   = new GammaMatrix[Nquark];
  GammaMatrix *sigma = new GammaMatrix[m_ndim * (m_ndim - 1) / 2 * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i] = gmset->get_GM(gmset->GAMMA5);
    gm0[i] = gmset->get_GM(gmset->UNITY);
    sigma[0 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA12);
    sigma[1 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA23);
    sigma[2 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA31);
    sigma[3 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA41);
    sigma[4 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA42);
    sigma[5 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA43);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  corr.meson_correlator_x(mcorr, gm0[0], gm0[0], cq, cq);
  for (int t = 0; t < Lx; ++t) {
    vout.general(m_vl, "CC_xcorr %lf %4d %0.16e %0.16e\n",
                 tt, t, real(mcorr[t]), imag(mcorr[t]));
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator_x(mcorr, gm0[iq], gm0[iq], sq[iq], cq);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "CS_xcorr_ss_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
    corr.meson_correlator_x(mcorr, gm0[iq], gm0[iq], cq, sq[iq]);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "SC_xcorr_ss_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator_x(mcorr, gm5[iq], gm5[iq], sq[iq], cq);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "CS_xcorr_pp_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
    corr.meson_correlator_x(mcorr, gm5[iq], gm5[iq], cq, sq[iq]);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "SC_xcorr_pp_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], sq[iq], cq);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "CS_xcorr_vv_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
      corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], cq, sq[iq]);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "SC_xcorr_vv_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator_x(mcorr, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], sq[iq], cq);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "CS_xcorr_aa_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
      corr.meson_correlator_x(mcorr, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], cq, sq[iq]);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "SC_xcorr_aa_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    corr.meson_correlator_x(mcorr, g5mu[m_ndim * iq], gm5[iq], sq[iq], cq);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "CS_xcorr_a0p_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
    corr.meson_correlator_x(mcorr, g5mu[m_ndim * iq], gm5[iq], cq, sq[iq]);
    for (int t = 0; t < Lx; ++t) {
      vout.general(m_vl, "SC_xcorr_a0p_iq%d %lf %4d %0.16e %0.16e\n",
                   iq, tt, t, real(mcorr[t]), imag(mcorr[t]));
    }
  }

  for (int mu = 0; mu < m_ndim * (m_ndim - 1) / 2; ++mu) {
    for (int iq = 0; iq < Nquark; ++iq) {
      corr.meson_correlator_x(mcorr, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sq[iq], cq);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "CS_xcorr_tt_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
      corr.meson_correlator_x(mcorr, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], cq, sq[iq]);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "SC_xcorr_tt_iq%d %d %lf %4d %0.16e %0.16e\n",
                     iq, mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] sigma;
  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_meson_correlator_t_FT(double tt,
                                                             const vector<vector<Field_F> >& sq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = sq.size();

  assert(Nquark == 2);

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);
  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm    = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *g5mu  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm5   = new GammaMatrix[Nquark];
  GammaMatrix *gm0   = new GammaMatrix[Nquark];
  GammaMatrix *sigma = new GammaMatrix[m_ndim * (m_ndim - 1) / 2 * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i] = gmset->get_GM(gmset->GAMMA5);
    gm0[i] = gmset->get_GM(gmset->UNITY);
    sigma[0 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA12);
    sigma[1 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA23);
    sigma[2 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA31);
    sigma[3 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA41);
    sigma[4 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA42);
    sigma[5 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA43);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator(mcorr, momentum_sink, gm0[iq], gm0[iq], sq[iq], sq[iq], source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "SS_tcorr_FT_ss_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
        corr.meson_momentum_correlator(mcorr, momentum_sink, gm0[1], gm0[0], sq[1], sq[0], source_position);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "SS_tcorr_FT_ss_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                       0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
        }
        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator(mcorr, momentum_sink, gm5[iq], gm5[iq], sq[iq], sq[iq], source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "SS_tcorr_FT_pp_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
        corr.meson_momentum_correlator(mcorr, momentum_sink, gm5[1], gm5[0], sq[1], sq[0], source_position);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "SS_tcorr_FT_pp_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                       0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
        }
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], sq[iq], sq[iq], source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "SS_tcorr_FT_vv_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
          corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim], gm[mu], sq[1], sq[0], source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "SS_tcorr_FT_vv_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator(mcorr, momentum_sink, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], sq[iq], sq[iq], source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "SS_tcorr_FT_aa_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
          corr.meson_momentum_correlator(mcorr, momentum_sink, g5mu[mu + m_ndim], g5mu[mu], sq[1], sq[0], source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "SS_tcorr_FT_aa_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator(mcorr, momentum_sink, g5mu[3 + m_ndim * iq], gm5[iq], sq[iq], sq[iq], source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "SS_tcorr_FT_a3p_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
        corr.meson_momentum_correlator(mcorr, momentum_sink, g5mu[3 + m_ndim], gm5[0], sq[1], sq[0], source_position);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "SS_tcorr_FT_a3p_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                       0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
        }
        for (int mu = 0; mu < m_ndim * (m_ndim - 1) / 2; ++mu) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator(mcorr, momentum_sink, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sq[iq], sq[iq], source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "SS_tcorr_FT_tt_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
          corr.meson_momentum_correlator(mcorr, momentum_sink, sigma[mu + m_ndim * (m_ndim - 1) / 2], sigma[mu], sq[1], sq[0], source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "SS_tcorr_FT_tt_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  delete [] sigma;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_meson_correlator_t_FT(double tt,
                                                             const vector<vector<Field_F> >& sq,
                                                             const vector<Field_F>& cq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = sq.size();

  assert(Nquark > 1);

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);
  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm    = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *g5mu  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm5   = new GammaMatrix[Nquark];
  GammaMatrix *gm0   = new GammaMatrix[Nquark];
  GammaMatrix *sigma = new GammaMatrix[m_ndim * (m_ndim - 1) / 2 * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i] = gmset->get_GM(gmset->GAMMA5);
    gm0[i] = gmset->get_GM(gmset->UNITY);
    sigma[0 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA12);
    sigma[1 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA23);
    sigma[2 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA31);
    sigma[3 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA41);
    sigma[4 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA42);
    sigma[5 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA43);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        corr.meson_momentum_correlator(mcorr, momentum_sink, gm0[0], gm0[0], cq, cq, source_position);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "CC_tcorr_FT %d %d %d %lf %4d %0.16e %0.16e\n",
                       momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
        }

        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator(mcorr, momentum_sink, gm0[iq], gm0[iq], sq[iq], cq, source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "CS_tcorr_FT_ss_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
          corr.meson_momentum_correlator(mcorr, momentum_sink, gm0[iq], gm0[iq], cq, sq[iq], source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "SC_tcorr_FT_ss_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }

        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator(mcorr, momentum_sink, gm5[iq], gm5[iq], sq[iq], cq, source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "CS_tcorr_FT_pp_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
          corr.meson_momentum_correlator(mcorr, momentum_sink, gm5[iq], gm5[iq], cq, sq[iq], source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "SC_tcorr_FT_pp_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }

        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], sq[iq], cq, source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "CS_tcorr_FT_vv_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
            corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], cq, sq[iq], source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "SC_tcorr_FT_vv_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator(mcorr, momentum_sink, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], sq[iq], cq, source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "CS_tcorr_FT_aa_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
            corr.meson_momentum_correlator(mcorr, momentum_sink, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], cq, sq[iq], source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "SC_tcorr_FT_aa_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
        }

        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator(mcorr, momentum_sink, g5mu[3 + m_ndim * iq], gm5[iq], sq[iq], cq, source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "CS_tcorr_FT_a3p_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
          corr.meson_momentum_correlator(mcorr, momentum_sink, g5mu[3 + m_ndim * iq], gm5[iq], cq, sq[iq], source_position);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "SC_tcorr_FT_a3p_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }

        for (int mu = 0; mu < m_ndim * (m_ndim - 1) / 2; ++mu) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator(mcorr, momentum_sink, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sq[iq], cq, source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "CS_tcorr_FT_tt_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
            corr.meson_momentum_correlator(mcorr, momentum_sink, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], cq, sq[iq], source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "SC_tcorr_FT_tt_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  delete [] sigma;
  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  delete gmset;

  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_meson_correlator_x_FT(double tt,
                                                             const vector<vector<Field_F> >& sq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = sq.size();

  assert(Nquark == 2);

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);
  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm    = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *g5mu  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm5   = new GammaMatrix[Nquark];
  GammaMatrix *gm0   = new GammaMatrix[Nquark];
  GammaMatrix *sigma = new GammaMatrix[m_ndim * (m_ndim - 1) / 2 * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i] = gmset->get_GM(gmset->GAMMA5);
    gm0[i] = gmset->get_GM(gmset->UNITY);
    sigma[0 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA12);
    sigma[1 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA23);
    sigma[2 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA31);
    sigma[3 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA41);
    sigma[4 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA42);
    sigma[5 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA43);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //mentum_sink[0]=1;
        //mentum_sink[1]=0;
        //mentum_sink[2]=0;
        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm0[iq], gm0[iq], sq[iq], sq[iq], source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "SS_xcorr_FT_ss_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
        corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm0[1], gm0[0], sq[1], sq[0], source_position);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "SS_xcorr_FT_ss_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                       0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
        }
        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm5[iq], gm5[iq], sq[iq], sq[iq], source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "SS_xcorr_FT_pp_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
        corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm5[1], gm5[0], sq[1], sq[0], source_position);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "SS_xcorr_FT_pp_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                       0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
        }
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], sq[iq], sq[iq], source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "SS_xcorr_FT_vv_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim], gm[mu], sq[1], sq[0], source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "SS_xcorr_FT_vv_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], sq[iq], sq[iq], source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "SS_xcorr_FT_aa_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, g5mu[mu + m_ndim], g5mu[mu], sq[1], sq[0], source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "SS_xcorr_FT_aa_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, g5mu[m_ndim * iq], gm5[iq], sq[iq], sq[iq], source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "SS_xcorr_FT_a0p_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
        corr.meson_momentum_correlator_x(mcorr, momentum_sink, g5mu[m_ndim], gm5[0], sq[1], sq[0], source_position);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "SS_xcorr_FT_a0p_iq%d%d %d %d %d %lf %4d %0.16e %0.16e\n",
                       0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
        }
        for (int mu = 0; mu < m_ndim * (m_ndim - 1) / 2; ++mu) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sq[iq], sq[iq], source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "SS_xcorr_FT_tt_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, sigma[mu + m_ndim * (m_ndim - 1) / 2], sigma[mu], sq[1], sq[0], source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "SS_xcorr_FT_tt_iq%d%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         0, 1, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  delete [] sigma;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_meson_correlator_x_FT(double tt,
                                                             const vector<vector<Field_F> >& sq,
                                                             const vector<Field_F>& cq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = sq.size();

  assert(Nquark > 1);

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);
  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm    = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *g5mu  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm5   = new GammaMatrix[Nquark];
  GammaMatrix *gm0   = new GammaMatrix[Nquark];
  GammaMatrix *sigma = new GammaMatrix[m_ndim * (m_ndim - 1) / 2 * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i] = gmset->get_GM(gmset->GAMMA5);
    gm0[i] = gmset->get_GM(gmset->UNITY);
    sigma[0 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA12);
    sigma[1 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA23);
    sigma[2 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA31);
    sigma[3 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA41);
    sigma[4 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA42);
    sigma[5 + i * m_ndim * (m_ndim - 1) / 2] = gmset->get_GM(gmset->SIGMA43);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm0[0], gm0[0], cq, cq, source_position);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "CC_xcorr_FT %d %d %d %lf %4d %0.16e %0.16e\n",
                       momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
        }

        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm0[iq], gm0[iq], sq[iq], cq, source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "CS_xcorr_FT_ss_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm0[iq], gm0[iq], cq, sq[iq], source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "SC_xcorr_FT_ss_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }

        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm5[iq], gm5[iq], sq[iq], cq, source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "CS_xcorr_FT_pp_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm5[iq], gm5[iq], cq, sq[iq], source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "SC_xcorr_FT_pp_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }

        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], sq[iq], cq, source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "CS_xcorr_FT_vv_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[mu + m_ndim * iq], cq, sq[iq], source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "SC_xcorr_FT_vv_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], sq[iq], cq, source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "CS_xcorr_FT_aa_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, g5mu[mu + m_ndim * iq], g5mu[mu + m_ndim * iq], cq, sq[iq], source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "SC_xcorr_FT_aa_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
        }

        for (int iq = 0; iq < Nquark; ++iq) {
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, g5mu[m_ndim * iq], gm5[iq], sq[iq], cq, source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "CS_xcorr_FT_a0p_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
          corr.meson_momentum_correlator_x(mcorr, momentum_sink, g5mu[m_ndim * iq], gm5[iq], cq, sq[iq], source_position);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "SC_xcorr_FT_a0p_iq%d %d %d %d %lf %4d %0.16e %0.16e\n",
                         iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, real(mcorr[t]), imag(mcorr[t]));
          }
        }

        for (int mu = 0; mu < m_ndim * (m_ndim - 1) / 2; ++mu) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sq[iq], cq, source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "CS_xcorr_FT_tt_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], sigma[mu + m_ndim * (m_ndim - 1) / 2 * iq], cq, sq[iq], source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "SC_xcorr_FT_tt_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, t, real(mcorr[t]), imag(mcorr[t]));
            }
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  delete [] sigma;
  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  delete gmset;

  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O3a_correlator_t(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = m_params_clover->size();

  GammaMatrix *gm = new GammaMatrix[m_ndim * Nquark];

  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int nu = 0; nu < m_ndim; ++nu) {
      for (int rho = 0; rho < m_ndim; ++rho) {
        for (int sigma = 0; sigma < m_ndim; ++sigma) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dsq[nu][0][iq], sq[sigma + 1][iq]);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "O3O3a_SS_tcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
            }
            corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dsq[nu][0][iq], cq[sigma + 1]);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "O3O3a_SC_tcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
            }
            corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dcq[nu][0], sq[sigma + 1][iq]);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "O3O3a_CS_tcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
            }
          }
          corr.meson_correlator(mcorr, gm[mu], gm[rho], dcq[nu][0], cq[sigma + 1]);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "O3O3a_CC_tcorr %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O3b_correlator_t(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = m_params_clover->size();

  GammaMatrix *gm = new GammaMatrix[m_ndim * Nquark];

  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int nu = 0; nu < m_ndim; ++nu) {
      for (int rho = 0; rho < m_ndim; ++rho) {
        for (int sigma = 0; sigma < m_ndim; ++sigma) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], sq[0][iq], dsq[nu][sigma + 1][iq]);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "O3O3b_SS_tcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
            }
            corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], sq[0][iq], dcq[nu][sigma + 1]);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "O3O3b_SC_tcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
            }
            corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], cq[0], dsq[nu][sigma + 1][iq]);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "O3O3b_CS_tcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
            }
          }
          corr.meson_correlator(mcorr, gm[mu], gm[rho], cq[0], dcq[nu][sigma + 1]);
          for (int t = 0; t < Lt; ++t) {
            vout.general(m_vl, "O3O3b_CC_tcorr %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O5_correlator_t(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = m_params_clover->size();

  GammaMatrix *gm  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm0 = new GammaMatrix[Nquark];

  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int nu = 0; nu < m_ndim; ++nu) {
      for (int iq = 0; iq < Nquark; ++iq) {
        corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm0[iq], dsq[nu][0][iq], sq[0][iq]);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "O3O5_SS_tcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
        }
        corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm0[iq], dsq[nu][0][iq], cq[0]);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "O3O5_SC_tcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
        }
        corr.meson_correlator(mcorr, gm[mu + m_ndim * iq], gm0[iq], dcq[nu][0], sq[0][iq]);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "O3O5_CS_tcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
        }
      }
      corr.meson_correlator(mcorr, gm[mu], gm0[0], dcq[nu][0], cq[0]);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "O3O5_CC_tcorr %d %d %lf %4d %0.16e %0.16e\n",
                     mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] gm0;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O5O3_correlator_t(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = m_params_clover->size();

  GammaMatrix *gm  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm0 = new GammaMatrix[Nquark];

  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int nu = 0; nu < m_ndim; ++nu) {
      for (int iq = 0; iq < Nquark; ++iq) {
        corr.meson_correlator(mcorr, gm0[iq], gm[mu + m_ndim * iq], sq[nu + 1][iq], sq[0][iq]);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "O5O3_SS_tcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
        }
        corr.meson_correlator(mcorr, gm0[iq], gm[mu + m_ndim * iq], sq[nu + 1][iq], cq[0]);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "O5O3_SC_tcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
        }
        corr.meson_correlator(mcorr, gm0[iq], gm[mu + m_ndim * iq], cq[nu + 1], sq[0][iq]);
        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "O5O3_CS_tcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
        }
      }
      corr.meson_correlator(mcorr, gm0[0], gm[mu], cq[nu + 1], cq[0]);
      for (int t = 0; t < Lt; ++t) {
        vout.general(m_vl, "O5O3_CC_tcorr %d %d %lf %4d %0.16e %0.16e\n",
                     mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] gm0;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O3a_correlator_x(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = m_params_clover->size();

  GammaMatrix *gm = new GammaMatrix[m_ndim * Nquark];

  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int nu = 0; nu < m_ndim; ++nu) {
      for (int rho = 0; rho < m_ndim; ++rho) {
        for (int sigma = 0; sigma < m_ndim; ++sigma) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dsq[nu][0][iq], sq[sigma + 1][iq]);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "O3O3a_SS_xcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
            }
            corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dsq[nu][0][iq], cq[sigma + 1]);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "O3O3a_SC_xcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
            }
            corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dcq[nu][0], sq[sigma + 1][iq]);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "O3O3a_CS_xcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
            }
          }
          corr.meson_correlator_x(mcorr, gm[mu], gm[rho], dcq[nu][0], cq[sigma + 1]);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "O3O3a_CC_xcorr %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O3b_correlator_x(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = m_params_clover->size();

  GammaMatrix *gm = new GammaMatrix[m_ndim * Nquark];

  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int nu = 0; nu < m_ndim; ++nu) {
      for (int rho = 0; rho < m_ndim; ++rho) {
        for (int sigma = 0; sigma < m_ndim; ++sigma) {
          for (int iq = 0; iq < Nquark; ++iq) {
            corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], sq[0][iq], dsq[nu][sigma + 1][iq]);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "O3O3b_SS_xcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
            }
            corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], sq[0][iq], dcq[nu][sigma + 1]);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "O3O3b_SC_xcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
            }
            corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], cq[0], dsq[nu][sigma + 1][iq]);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "O3O3b_CS_xcorr_iq%d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           iq, mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
            }
          }
          corr.meson_correlator_x(mcorr, gm[mu], gm[rho], cq[0], dcq[nu][sigma + 1]);
          for (int t = 0; t < Lx; ++t) {
            vout.general(m_vl, "O3O3b_CC_xcorr %d %d %d %d %lf %4d %0.16e %0.16e\n",
                         mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O5_correlator_x(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = m_params_clover->size();

  GammaMatrix *gm  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm0 = new GammaMatrix[Nquark];

  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int nu = 0; nu < m_ndim; ++nu) {
      for (int iq = 0; iq < Nquark; ++iq) {
        corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm0[iq], dsq[nu][0][iq], sq[0][iq]);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "O3O5_SS_xcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
        }
        corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm0[iq], dsq[nu][0][iq], cq[0]);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "O3O5_SC_xcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
        }
        corr.meson_correlator_x(mcorr, gm[mu + m_ndim * iq], gm0[iq], dcq[nu][0], sq[0][iq]);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "O3O5_CS_xcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
        }
      }
      corr.meson_correlator_x(mcorr, gm[mu], gm0[0], dcq[nu][0], cq[0]);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "O3O5_CC_xcorr %d %d %lf %4d %0.16e %0.16e\n",
                     mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
      }
    }
  }
  double result = real(mcorr[0]);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  delete [] gm;
  delete [] gm0;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O5O3_correlator_x(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = m_params_clover->size();

  GammaMatrix *gm  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm0 = new GammaMatrix[Nquark];

  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < m_ndim; ++mu) {
    for (int nu = 0; nu < m_ndim; ++nu) {
      for (int iq = 0; iq < Nquark; ++iq) {
        corr.meson_correlator_x(mcorr, gm0[iq], gm[mu + m_ndim * iq], sq[nu + 1][iq], sq[0][iq]);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "O5O3_SS_xcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
        }
        corr.meson_correlator_x(mcorr, gm0[iq], gm[mu + m_ndim * iq], sq[nu + 1][iq], cq[0]);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "O5O3_SC_xcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
        }
        corr.meson_correlator_x(mcorr, gm0[iq], gm[mu + m_ndim * iq], cq[nu + 1], sq[0][iq]);
        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "O5O3_CS_xcorr_iq%d %d %d %lf %4d %0.16e %0.16e\n",
                       iq, mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
        }
      }
      corr.meson_correlator_x(mcorr, gm0[0], gm[mu], cq[nu + 1], cq[0]);
      for (int t = 0; t < Lx; ++t) {
        vout.general(m_vl, "O5O3_CC_xcorr %d %d %lf %4d %0.16e %0.16e\n",
                     mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
      }
    }
  }
  double result = real(mcorr[0]);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  delete [] gm;
  delete [] gm0;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O3a_correlator_t_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = m_params_clover->size();

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);

  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm = new GammaMatrix[m_ndim * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int nu = 0; nu < m_ndim; ++nu) {
            for (int rho = 0; rho < m_ndim; ++rho) {
              for (int sigma = 0; sigma < m_ndim; ++sigma) {
                for (int iq = 0; iq < Nquark; ++iq) {
                  corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dsq[nu][0][iq], sq[sigma + 1][iq], source_position);
                  for (int t = 0; t < Lt; ++t) {
                    vout.general(m_vl, "O3O3a_SS_tcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
                  }
                  corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dsq[nu][0][iq], cq[sigma + 1], source_position);
                  for (int t = 0; t < Lt; ++t) {
                    vout.general(m_vl, "O3O3a_SC_tcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
                  }
                  corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dcq[nu][0], sq[sigma + 1][iq], source_position);
                  for (int t = 0; t < Lt; ++t) {
                    vout.general(m_vl, "O3O3a_CS_tcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
                  }
                }
                corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu], gm[rho], dcq[nu][0], cq[sigma + 1], source_position);
                for (int t = 0; t < Lt; ++t) {
                  vout.general(m_vl, "O3O3a_CC_tcorr_FT %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                               momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
                }
              }
            }
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O3b_correlator_t_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = m_params_clover->size();

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);

  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm = new GammaMatrix[m_ndim * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int nu = 0; nu < m_ndim; ++nu) {
            for (int rho = 0; rho < m_ndim; ++rho) {
              for (int sigma = 0; sigma < m_ndim; ++sigma) {
                for (int iq = 0; iq < Nquark; ++iq) {
                  corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], sq[0][iq], dsq[nu][sigma + 1][iq], source_position);
                  for (int t = 0; t < Lt; ++t) {
                    vout.general(m_vl, "O3O3b_SS_tcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
                  }
                  corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], sq[0][iq], dcq[nu][sigma + 1], source_position);
                  for (int t = 0; t < Lt; ++t) {
                    vout.general(m_vl, "O3O3b_SC_tcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
                  }
                  corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], cq[0], dsq[nu][sigma + 1][iq], source_position);
                  for (int t = 0; t < Lt; ++t) {
                    vout.general(m_vl, "O3O3b_CS_tcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
                  }
                }
                corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu], gm[rho], cq[0], dcq[nu][sigma + 1], source_position);
                for (int t = 0; t < Lt; ++t) {
                  vout.general(m_vl, "O3O3b_CC_tcorr_FT %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                               momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
                }
              }
            }
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O5_correlator_t_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = m_params_clover->size();

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);

  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm0 = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int nu = 0; nu < m_ndim; ++nu) {
            for (int iq = 0; iq < Nquark; ++iq) {
              corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm0[iq], dsq[nu][0][iq], sq[0][iq], source_position);
              for (int t = 0; t < Lt; ++t) {
                vout.general(m_vl, "O3O5_SS_tcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
              }
              corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm0[iq], dsq[nu][0][iq], cq[0], source_position);
              for (int t = 0; t < Lt; ++t) {
                vout.general(m_vl, "O3O5_SC_tcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
              }
              corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm0[iq], dcq[nu][0], sq[0][iq], source_position);
              for (int t = 0; t < Lt; ++t) {
                vout.general(m_vl, "O3O5_CS_tcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
              }
            }
            corr.meson_momentum_correlator(mcorr, momentum_sink, gm[mu], gm0[0], dcq[nu][0], cq[0], source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "O3O5_CC_tcorr_FT %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
            }
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] gm0;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O5O3_correlator_t_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq)
{
  int Lt     = CommonParameters::Lt();
  int Nquark = m_params_clover->size();

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);

  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm0 = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lt);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int nu = 0; nu < m_ndim; ++nu) {
            for (int iq = 0; iq < Nquark; ++iq) {
              corr.meson_momentum_correlator(mcorr, momentum_sink, gm0[iq], gm[mu + m_ndim * iq], sq[nu + 1][iq], sq[0][iq], source_position);
              for (int t = 0; t < Lt; ++t) {
                vout.general(m_vl, "O5O3_SS_tcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
              }
              corr.meson_momentum_correlator(mcorr, momentum_sink, gm0[iq], gm[mu + m_ndim * iq], sq[nu + 1][iq], cq[0], source_position);
              for (int t = 0; t < Lt; ++t) {
                vout.general(m_vl, "O5O3_SC_tcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
              }
              corr.meson_momentum_correlator(mcorr, momentum_sink, gm0[iq], gm[mu + m_ndim * iq], cq[nu + 1], sq[0][iq], source_position);
              for (int t = 0; t < Lt; ++t) {
                vout.general(m_vl, "O5O3_CS_tcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
              }
            }
            corr.meson_momentum_correlator(mcorr, momentum_sink, gm0[0], gm[mu], cq[nu + 1], cq[0], source_position);
            for (int t = 0; t < Lt; ++t) {
              vout.general(m_vl, "O5O3_CC_tcorr_FT %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
            }
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] gm0;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O3a_correlator_x_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = m_params_clover->size();

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);

  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm = new GammaMatrix[m_ndim * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int nu = 0; nu < m_ndim; ++nu) {
            for (int rho = 0; rho < m_ndim; ++rho) {
              for (int sigma = 0; sigma < m_ndim; ++sigma) {
                for (int iq = 0; iq < Nquark; ++iq) {
                  corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dsq[nu][0][iq], sq[sigma + 1][iq], source_position);
                  for (int t = 0; t < Lx; ++t) {
                    vout.general(m_vl, "O3O3a_SS_xcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
                  }
                  corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dsq[nu][0][iq], cq[sigma + 1], source_position);
                  for (int t = 0; t < Lx; ++t) {
                    vout.general(m_vl, "O3O3a_SC_xcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
                  }
                  corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], dcq[nu][0], sq[sigma + 1][iq], source_position);
                  for (int t = 0; t < Lx; ++t) {
                    vout.general(m_vl, "O3O3a_CS_xcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
                  }
                }
                corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu], gm[rho], dcq[nu][0], cq[sigma + 1], source_position);
                for (int t = 0; t < Lx; ++t) {
                  vout.general(m_vl, "O3O3a_CC_xcorr_FT %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                               momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
                }
              }
            }
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O3b_correlator_x_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = m_params_clover->size();

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);

  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm = new GammaMatrix[m_ndim * Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int nu = 0; nu < m_ndim; ++nu) {
            for (int rho = 0; rho < m_ndim; ++rho) {
              for (int sigma = 0; sigma < m_ndim; ++sigma) {
                for (int iq = 0; iq < Nquark; ++iq) {
                  corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], sq[0][iq], dsq[nu][sigma + 1][iq], source_position);
                  for (int t = 0; t < Lx; ++t) {
                    vout.general(m_vl, "O3O3b_SS_xcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
                  }
                  corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], sq[0][iq], dcq[nu][sigma + 1], source_position);
                  for (int t = 0; t < Lx; ++t) {
                    vout.general(m_vl, "O3O3b_SC_xcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
                  }
                  corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm[rho + m_ndim * iq], cq[0], dsq[nu][sigma + 1][iq], source_position);
                  for (int t = 0; t < Lx; ++t) {
                    vout.general(m_vl, "O3O3b_CS_xcorr_FT_iq%d %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                                 iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
                  }
                }
                corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu], gm[rho], cq[0], dcq[nu][sigma + 1], source_position);
                for (int t = 0; t < Lx; ++t) {
                  vout.general(m_vl, "O3O3b_CC_xcorr_FT %d %d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                               momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, rho, sigma, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
                }
              }
            }
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O3O5_correlator_x_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq, const std::vector<std::vector<std::vector<std::vector<Field_F> > > >& dsq, const std::vector<std::vector<std::vector<Field_F> > >& dcq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = m_params_clover->size();

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);

  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm0 = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int nu = 0; nu < m_ndim; ++nu) {
            for (int iq = 0; iq < Nquark; ++iq) {
              corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm0[iq], dsq[nu][0][iq], sq[0][iq], source_position);
              for (int t = 0; t < Lx; ++t) {
                vout.general(m_vl, "O3O5_SS_xcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
              }
              corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm0[iq], dsq[nu][0][iq], cq[0], source_position);
              for (int t = 0; t < Lx; ++t) {
                vout.general(m_vl, "O3O5_SC_xcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
              }
              corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu + m_ndim * iq], gm0[iq], dcq[nu][0], sq[0][iq], source_position);
              for (int t = 0; t < Lx; ++t) {
                vout.general(m_vl, "O3O5_CS_xcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
              }
            }
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm[mu], gm0[0], dcq[nu][0], cq[0], source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "O3O5_CC_xcorr_FT %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, -2 * real(mcorr[t]), -2 * imag(mcorr[t]));
            }
          }
        }
        if (m_filename_output != "stdout") {
          log_file.close();
          vout.init(log_file_previous);
        }
      }
    }
  }
  double result = real(mcorr[0]);
  delete [] gm;
  delete [] gm0;
  delete gmset;
  return result;
}


//====================================================================
double FermionFlow_2pt_Function::print_O5O3_correlator_x_FT(double tt, const std::vector<std::vector<std::vector<Field_F> > >& sq, const std::vector<std::vector<Field_F> >& cq)
{
  int Lx     = CommonParameters::Lx();
  int Nquark = m_params_clover->size();

  //m_max_mom=0;
  int Np = (2 * m_max_mom + 1);

  // vector<int> source_position(3, 0);
  // vector<int> momentum_sink(3);
  vector<int> source_position(m_ndim, 0);
  vector<int> momentum_sink(m_ndim - 1);

  GammaMatrix *gm  = new GammaMatrix[m_ndim * Nquark];
  GammaMatrix *gm0 = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * m_ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * m_ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * m_ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * m_ndim] = gmset->get_GM(gmset->GAMMA4);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }
  m_str_gmset_type = (*m_params_clover)[0].get_string("gamma_matrix_type");
  GammaMatrixSet   *gmset = GammaMatrixSet::New(m_str_gmset_type);
  Corr2pt_4spinor  corr(gmset);
  vector<dcomplex> mcorr(Lx);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        //momentum_sink[0]=1;
        //momentum_sink[1]=0;
        //momentum_sink[2]=0;
        for (int mu = 0; mu < m_ndim; ++mu) {
          for (int nu = 0; nu < m_ndim; ++nu) {
            for (int iq = 0; iq < Nquark; ++iq) {
              corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm0[iq], gm[mu + m_ndim * iq], sq[nu + 1][iq], sq[0][iq], source_position);
              for (int t = 0; t < Lx; ++t) {
                vout.general(m_vl, "O5O3_SS_xcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
              }
              corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm0[iq], gm[mu + m_ndim * iq], sq[nu + 1][iq], cq[0], source_position);
              for (int t = 0; t < Lx; ++t) {
                vout.general(m_vl, "O5O3_SC_xcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
              }
              corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm0[iq], gm[mu + m_ndim * iq], cq[nu + 1], sq[0][iq], source_position);
              for (int t = 0; t < Lx; ++t) {
                vout.general(m_vl, "O5O3_CS_xcorr_FT_iq%d %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                             iq, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
              }
            }
            corr.meson_momentum_correlator_x(mcorr, momentum_sink, gm0[0], gm[mu], cq[nu + 1], cq[0], source_position);
            for (int t = 0; t < Lx; ++t) {
              vout.general(m_vl, "O5O3_CC_xcorr_FT %d %d %d %d %d %lf %4d %0.16e %0.16e\n",
                           momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, t, 2 * real(mcorr[t]), 2 * imag(mcorr[t]));
            }
          }
        }
      }
    }
  }
  double result = real(mcorr[0]);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] gm0;
  delete gmset;
  return result;
}


//====================================================================
//============================================================END=====
