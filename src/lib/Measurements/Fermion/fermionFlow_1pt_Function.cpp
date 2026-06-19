/*!
        @file    fermionFlow_1pt_Function.cpp

        @brief

        @author  <Yusuke Taniguchi> tanigchi@het.ph.tsukuba.ac.jp
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#include "fermionFlow_1pt_Function.h"

const std::string FermionFlow_1pt_Function::class_name = "FermionFlow_1pt_Function";

//====================================================================
void FermionFlow_1pt_Function::set_parameters(const Parameters& params_measurement, const Parameters& params_gflow, const Parameters& params_fflow, const Parameters& params_source_random)
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
  string noise_type;

  int err = 0;
  err += params_measurement.fetch_int("number_of_noises", Nnoise);
  err += params_measurement.fetch_int("initial_tau", initial_tau);
  err += params_measurement.fetch_int("number_of_measurement_times", number_of_measurement_times);
  err += params_measurement.fetch_int("measurement_interval", measurement_interval);
  err += params_measurement.fetch_int("gauge_store_interval", gauge_store_interval);
  err += params_measurement.fetch_int("max_momentum", max_mom);
  err += params_fflow.fetch_double("step_size", step_size);
  err += params_source_random.fetch_string("noise_type", noise_type);

  if (err) {
    vout.crucial(m_vl, "%s: fetch error, input parameter not found.\n", class_name.c_str());
    abort();
  }

  assert(number_of_measurement_times >= 1);
  assert(initial_tau > 0);

  set_parameters(Nnoise, number_of_measurement_times, measurement_interval, gauge_store_interval, initial_tau, step_size, max_mom);
  m_gflow.set_parameters(params_gflow);
  m_fflow.set_parameters(params_fflow);
  m_nv.set_parameters(params_source_random);
}


//====================================================================
void FermionFlow_1pt_Function::get_parameters(Parameters& params) const
{
  params.set_int("number_of_noises", m_Nnoise);
  params.set_int("initial_tau", m_initial_tau);
  params.set_int("number_of_measurement_times", m_number_of_measurement_times);
  params.set_int("measurement_interval", m_measurement_interval);
  params.set_int("gauge_store_interval", m_gauge_store_interval);
  params.set_int("max_momentum", m_max_mom);
  params.set_double("step_size", m_step_size);
  // params.set_string("noise_type", m_noise_type);

  Parameters params_gflow;
  m_gflow.get_parameters(params_gflow);
  params.set_Parameters("gradient_flow", params_gflow);

  Parameters params_fflow;
  m_fflow.get_parameters(params_fflow);
  params.set_Parameters("fermion_flow", params_fflow);

  params.set_string("filename_output", m_filename_output);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  return;
}


//====================================================================
void FermionFlow_1pt_Function::set_parameters(const int Nnoise, const int number_of_measurement_times, const int measurement_interval, const int gauge_store_interval, const int initial_tau, const double step_size, const int max_mom)
{
  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  Nnoise = %d\n", Nnoise);
  vout.general(m_vl, "  initial_tau= %d\n", initial_tau);
  vout.general(m_vl, "  number_of_measurement_times = %d\n", number_of_measurement_times);
  vout.general(m_vl, "  measurement_interval = %d\n", measurement_interval);
  vout.general(m_vl, "  gauge_store_interval = %d\n", gauge_store_interval);
  vout.general(m_vl, "  step_size = %lf\n", step_size);
  vout.general(m_vl, "  max_momentum = %d\n", max_mom);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Nnoise);

  err += ParameterCheck::is_satisfied(number_of_measurement_times >= 1);
  err += ParameterCheck::is_satisfied(initial_tau > 0);
  // assert(number_of_measurement_times>=1);
  // assert(initial_tau>0);

  if (err) {
    vout.crucial(m_vl, "%s: parameter range check failed.\n", class_name.c_str());
    abort();
  }

  //- store values
  m_Nnoise = Nnoise;
  m_number_of_measurement_times = number_of_measurement_times;
  m_measurement_interval        = measurement_interval;
  m_gauge_store_interval        = gauge_store_interval;
  m_step_size   = step_size;
  m_initial_tau = initial_tau;
  m_max_mom     = max_mom;
}


//====================================================================
double FermionFlow_1pt_Function::measure_adjoint_flow(Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  assert(m_gauge_store_interval > 0);
  assert(m_number_of_measurement_times >= m_gauge_store_interval);
  assert(m_number_of_measurement_times % m_gauge_store_interval == 0);
  assert(m_number_of_measurement_times >= 1);
  assert(m_initial_tau > 0);

  double t; //dummy
  //Number of confs at non-zero flow time to store in U0
  int number_of_confs = (m_initial_tau + m_number_of_measurement_times - 1) / m_gauge_store_interval;

  Field_G U0(Nvol, Ndim *(number_of_confs + 1));

  // Store Gauge field at flow time tau=0
  for (int mu = 0; mu < Ndim; ++mu) {
    U0.setpart_ex(mu, U, mu);
  }
  for (int i = 0; i < number_of_confs; ++i) {
    // Store Gauge field at flow time tau
    for (int j = 0; j < m_gauge_store_interval * m_measurement_interval; j++) {
      m_gflow.evolve(t, U);
    }
    for (int mu = 0; mu < Ndim; ++mu) {
      U0.setpart_ex(mu + Ndim * (i + 1), U, mu);
    }
  }

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Number of noise vector = %d\n", m_Nnoise);

  int    Nconv;
  double diff;

  const int Nquark  = m_fprop_lex.size();
  const int Nquark2 = m_params_clover->size();
  assert(Nquark == Nquark2);

  const int Nex = 1;
  const int Nin = 2 * m_nc * m_nd;

  Field   xi(Nin, Nvol, Nex);
  Field   v(Nin, Nvol, Nex);
  Field_F w;
  Field_F psi;

  Field_F_1spinor eta(Nvol, Nex);
  Field_F_1spinor v2(Nvol, Nex);

  Field_F_1spinor xi2(Nvol, (1 + Ndim) * m_nc * m_Nnoise * Nex);

  vector<dcomplex> tr_scalar;
  vector<dcomplex> tr_xidxi;
  vector<dcomplex> tr_tensor;
  vector<dcomplex> tr_pseudo;
  vector<dcomplex> tr_vector;
  vector<dcomplex> tr_axial;
  tr_scalar.resize(Nquark);
  tr_xidxi.resize(Nquark);
  tr_tensor.resize(Ndim * Ndim * Nquark);
  tr_pseudo.resize(Nquark);
  tr_vector.resize(Nquark * Ndim);
  tr_axial.resize(Nquark * Ndim);

  GammaMatrix *gm  = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *gm5 = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * Ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * Ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * Ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * Ndim] = gmset->get_GM(gmset->GAMMA4);
    gm5[i]           = gmset->get_GM(gmset->GAMMA5);
    delete gmset;
  }

  for (int measure = m_initial_tau; measure < m_initial_tau + m_number_of_measurement_times; ++measure) {
    double tt = m_step_size * measure * m_measurement_interval;

    for (int i = 0; i < Nquark; ++i) {
      tr_scalar[i] = 0;
      tr_xidxi[i]  = 0;
      tr_pseudo[i] = 0;
    }
    for (int i = 0; i < Ndim * Ndim * Nquark; ++i) {
      tr_tensor[i] = 0;
    }
    for (int i = 0; i < Ndim * Nquark; ++i) {
      tr_vector[i] = 0;
      tr_axial[i]  = 0;
    }

    for (int mu = 0; mu < Ndim; ++mu) {
      U.setpart_ex(mu, U0, mu + Ndim * (measure / m_gauge_store_interval));
    }
    if (measure % m_gauge_store_interval != 0) {
      for (int i = 0; i < (measure % m_gauge_store_interval) * m_measurement_interval; i++) {
        m_gflow.evolve(t, U);
      }
    }
    for (int inoise = 0; inoise < m_Nnoise; ++inoise) {
      for (int ic = 0; ic < m_nc; ++ic) {
        // noise vector was set.
        m_nv.set_all_space_time(eta, ic);
        for (int isite = 0; isite < xi2.nvol(); isite++) {
          for (int icc = 0; icc < m_nc; icc++) {
            int lex = iex(0, ic, inoise);
            xi2.set_ri(icc, isite, lex, eta.cmp_r(icc, isite, 0), eta.cmp_i(icc, isite, 0));
          }
        }
        for (int mu = 0; mu < Ndim; ++mu) {
          m_fflow.del_symmetric(mu, v2, eta, U);
          for (int isite = 0; isite < xi2.nvol(); isite++) {
            for (int icc = 0; icc < m_nc; icc++) {
              int lex = iex(mu + 1, ic, inoise);
              xi2.set_ri(icc, isite, lex, v2.cmp_r(icc, isite, 0), v2.cmp_i(icc, isite, 0));
            }
          }
        }
      }
    }

    m_fflow.set_parameters(measure * m_measurement_interval);
    m_fflow.evolve2(xi2, U0, measure, m_gauge_store_interval);

    for (int inoise = 0; inoise < m_Nnoise; ++inoise) {
      for (int id = 0; id < m_nd; id++) {
        for (int ic = 0; ic < m_nc; ++ic) {
          xi.set(0.0);
          for (int isite = 0; isite < xi2.nvol(); isite++) {
            for (int icc = 0; icc < m_nc; icc++) {
              int lex = iex(0, ic, inoise);
              xi.set(2 * (icc + m_nc * id), isite, 0, xi2.cmp_r(icc, isite, lex));
              xi.set(1 + 2 * (icc + m_nc * id), isite, 0, xi2.cmp_i(icc, isite, lex));
            }
          }
          for (int mu = 0; mu < Ndim; ++mu) {
            U.setpart_ex(mu, U0, mu);
          }
          for (int iq = 0; iq < Nquark; ++iq) {
            m_fprop_lex[iq]->invert_D(v, xi, Nconv, diff);
            vout.detailed(m_vl, "    Nconv = %d  diff  = %.8e\n", Nconv, diff);
            // now v is  M^-1 * xi.

            dcomplex tr = dotc(xi, v);
            tr_scalar[iq] += tr;
            vout.general(m_vl, " dxiSxi_flavor%d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, ic, id, tt, real(tr), imag(tr));
            tr            = dotc(xi, xi);
            tr_xidxi[iq] += tr;
            vout.general(m_vl, " dxixi_flavor%d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, ic, id, tt, real(tr), imag(tr));

            mult_GM(w, gm5[iq], xi);
            tr             = dotc(w, v);
            tr_pseudo[iq] += tr;
            vout.general(m_vl, " dxig5Sxi_flavor%d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, ic, id, tt, real(tr), imag(tr));

            for (int mu = 0; mu < Ndim; ++mu) {
              int k = mu + Ndim * iq;
              mult_GM(w, gm[mu + Ndim * iq], xi);
              tr            = dotc(w, v);
              tr_vector[k] += tr;
              vout.general(m_vl, " dxigmuSxi_flavor%d %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, ic, id, mu, tt, real(tr), imag(tr));
              mult_GM(psi, gm5[iq], w);
              tr           = dotc(psi, v);
              tr_axial[k] += tr;
              vout.general(m_vl, " dxigmug5Sxi_flavor%d %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, ic, id, mu, tt, real(tr), imag(tr));
            }

            for (int nu = 0; nu < Ndim; ++nu) {
              w.set(0.0);
              for (int isite = 0; isite < xi2.nvol(); isite++) {
                for (int icc = 0; icc < m_nc; icc++) {
                  int lex = iex(nu + 1, ic, inoise);
                  w.set_ri(icc, id, isite, 0, xi2.cmp_r(icc, isite, lex), xi2.cmp_i(icc, isite, lex));
                }
              }
              for (int mu = 0; mu < Ndim; ++mu) {
                int k = mu + Ndim * (nu + Ndim * iq);
                mult_GM(psi, gm[mu + Ndim * iq], w);
                tr            = dotc(psi, v);
                tr_tensor[k] += tr;
                vout.general(m_vl, " dpsiSxi_flavor%d %d %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, ic, id, mu, nu, tt, real(tr), imag(tr));
              }
            }
          }
        }
      }
    }
    for (int iq = 0; iq < Nquark; ++iq) {
      m_hopping_parameter = (*m_params_clover)[iq].get_double("hopping_parameter");
      tr_scalar[iq]      /= double(m_Nnoise) * Nvol * NPE;
      tr_scalar[iq]      *= 2.0 * m_hopping_parameter;
      tr_xidxi[iq]       /= double(m_Nnoise) * Nvol * NPE;
      vout.general(m_vl, " tr_scalar%d= %lf %0.16e %0.16e\n", iq, tt, -real(tr_scalar[iq]), -imag(tr_scalar[iq])); // (2.15) of uds.pdf
      vout.general(m_vl, " tr_xidxi%d= %lf %0.16e %0.16e\n", iq, tt, real(tr_xidxi[iq]), imag(tr_xidxi[iq]));      // (2.15) of uds.pdf

      tr_pseudo[iq] /= double(m_Nnoise) * Nvol * NPE;
      tr_pseudo[iq] *= 2.0 * m_hopping_parameter;
      vout.general(m_vl, " tr_pseudo%d= %lf %0.16e %0.16e\n", iq, tt, -real(tr_pseudo[iq]), -imag(tr_pseudo[iq])); // (2.15) of uds.pdf

      for (int nu = 0; nu < Ndim; ++nu) {
        for (int mu = 0; mu < Ndim; ++mu) {
          int k = mu + Ndim * (nu + Ndim * iq);
          tr_tensor[k] /= double(m_Nnoise) * Nvol * NPE;
          tr_tensor[k] *= 2.0 * m_hopping_parameter;
          vout.general(m_vl, " tr_tensor%d= %d %d %lf %0.16e %0.16e\n", iq, mu, nu, tt,
                       2 * real(tr_tensor[k]), 2 * imag(tr_tensor[k])); // (2.14) of uds.pdf
        }
      }

      for (int mu = 0; mu < Ndim; ++mu) {
        int k = mu + Ndim * iq;
        tr_vector[k] /= double(m_Nnoise) * Nvol * NPE;
        tr_vector[k] *= 2.0 * m_hopping_parameter;
        vout.general(m_vl, " tr_vector%d= %d %lf %0.16e %0.16e\n", iq, mu, tt, -real(tr_vector[k]), -imag(tr_vector[k]));
        tr_axial[k] /= double(m_Nnoise) * Nvol * NPE;
        tr_axial[k] *= 2.0 * m_hopping_parameter;
        vout.general(m_vl, " tr_axial%d= %d %lf %0.16e %0.16e\n", iq, mu, tt, -real(tr_axial[k]), -imag(tr_axial[k]));
      }
    }
  }
  delete [] gm;
  delete [] gm5;

  return 2 * real(tr_tensor[0]);
}


//====================================================================
double FermionFlow_1pt_Function::measure_normal_order(Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  double result = 0.0;

  assert(m_number_of_measurement_times >= 1);
  assert(m_initial_tau > 0);

  // Nstep=m_number_of_measurement_times*m_measurement_interval
  m_fflow.set_parameters(m_measurement_interval);

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Number of noise vector = %d\n", m_Nnoise);

  int    Nconv;
  double diff;

  Field_G U0((Field&)U);

  const int Nquark  = m_fprop_lex.size();
  const int Nquark2 = m_params_clover->size();
  assert(Nquark == Nquark2);

  const int Nex1 = m_nc * m_nd;
  const int Nex2 = m_nc;
  const int Nin  = 2 * m_nc * m_nd;

  Field_F         xi(Nvol, Nquark *Nex1);
  Field_F_1spinor xi2(Nvol, Nex2);

  Field v(Nin, Nvol, 1);
  Field w(Nin, Nvol, 1);

  Field_F x1;
  Field_F z1;

  Field_F_1spinor eta(Nvol, 1);

  /*
  GammaMatrix *gm = new GammaMatrix[Ndim*Nquark];
  GammaMatrix *gm5 = new GammaMatrix[Nquark];
  for(int i=0; i<Nquark; ++i){
    m_str_gmset_type=(*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0+i*Ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1+i*Ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2+i*Ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3+i*Ndim] = gmset->get_GM(gmset->GAMMA4);
    gm5[i] = gmset->get_GM(gmset->GAMMA5);
    delete gmset;
  }
  */
  for (int inoise = 0; inoise < m_Nnoise; ++inoise) {
    U = U0;
    for (int ic = 0; ic < m_nc; ++ic) {
      // set noide vector.
      m_nv.set_all_space_time(eta, ic);
      int lex = ic;
      for (int isite = 0; isite < Nvol; isite++) {
        for (int icc = 0; icc < m_nc; ++icc) {
          xi2.set_ri(icc, isite, lex, eta.cmp_r(icc, isite, 0), eta.cmp_i(icc, isite, 0));
        }
      }
      for (int id = 0; id < m_nd; ++id) {
        v.set(0.0);
        for (int isite = 0; isite < Nvol; isite++) {
          v.set(2 * (ic + m_nc * id), isite, 0, eta.cmp_r(ic, isite, 0));
          v.set(1 + 2 * (ic + m_nc * id), isite, 0, eta.cmp_i(ic, isite, 0));
        }
        for (int iq = 0; iq < Nquark; ++iq) {
          // w = M^-1 * v
          m_fprop_lex[iq]->invert_D(w, v, Nconv, diff);
          vout.general(m_vl, "   Nconv = %d  diff  = %.8e\n", Nconv, diff);
          lex = iq + Nquark * (ic + m_nc * id);
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
    for (int measure = m_initial_tau; measure < m_initial_tau + m_number_of_measurement_times; ++measure) {
      double tt = m_step_size * measure * m_measurement_interval;
      m_fflow.evolve_normal_order(xi2, xi, U);
      result += print_vev(inoise, tt, xi2, xi, U);
    }
  } // inoise
    //  delete [] gm;
    //  delete [] gm5;
  return result;
}


//====================================================================
double FermionFlow_1pt_Function::measure_disconnected(Field_G& U)
{
  const int Nvol   = CommonParameters::Nvol();
  double    result = 0.0;

  // Nstep=m_number_of_measurement_times*m_measurement_interval
  m_fflow.set_parameters(m_measurement_interval);

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Number of noise vector = %d\n", m_Nnoise);

  int    Nconv;
  double diff;

  Field_G U0((Field&)U);

  const int Nquark  = m_fprop_lex.size();
  const int Nquark2 = m_params_clover->size();
  assert(Nquark == Nquark2);

  const int Nex1 = m_nc * m_nd;
  const int Nex2 = m_nc;
  const int Nin  = 2 * m_nc * m_nd;

  Field_F         xi(Nvol, Nquark *Nex1);
  Field_F_1spinor xi2(Nvol, Nex2);

  Field v(Nin, Nvol, 1);
  Field w(Nin, Nvol, 1);

  Field_F_1spinor eta(Nvol, 1);


  for (int inoise = 0; inoise < m_Nnoise; ++inoise) {
    U = U0;
    for (int ic = 0; ic < m_nc; ++ic) {
      // set noide vector.
      m_nv.set_all_space_time(eta, ic);
      int lex = ic;
      for (int isite = 0; isite < Nvol; isite++) {
        for (int icc = 0; icc < m_nc; ++icc) {
          xi2.set_ri(icc, isite, lex, eta.cmp_r(icc, isite, 0), eta.cmp_i(icc, isite, 0));
        }
      }
      for (int id = 0; id < m_nd; ++id) {
        v.set(0.0);
        for (int isite = 0; isite < Nvol; isite++) {
          v.set(2 * (ic + m_nc * id), isite, 0, eta.cmp_r(ic, isite, 0));
          v.set(1 + 2 * (ic + m_nc * id), isite, 0, eta.cmp_i(ic, isite, 0));
        }
        for (int iq = 0; iq < Nquark; ++iq) {
          // w = M^-1 * v
          m_fprop_lex[iq]->invert_D(w, v, Nconv, diff);
          vout.general(m_vl, "   Nconv = %d  diff  = %.8e\n", Nconv, diff);
          lex = iq + Nquark * (ic + m_nc * id);
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
    for (int measure = m_initial_tau; measure < m_initial_tau + m_number_of_measurement_times; ++measure) {
      double tt = m_step_size * measure * m_measurement_interval;
      m_fflow.evolve_normal_order(xi2, xi, U);

      result += print_vev(inoise, tt, xi2, xi, U);
      result += print_correlator_t(inoise, tt, xi2, xi, U);
      result += print_correlator_t_FT(inoise, tt, xi2, xi, U);
      result += print_correlator_x(inoise, tt, xi2, xi, U);
      result += print_correlator_x_FT(inoise, tt, xi2, xi, U);

      /*
       print_correlator_y(inoise, tt, xi2, xi, U);
      print_correlator_y_FT(inoise, tt, xi2, xi, U);
      print_correlator_z(inoise, tt, xi2, xi, U);
      print_correlator_z_FT(inoise, tt, xi2, xi, U);
      */
    } // measure
  }   // inoise
  return result;
}


//====================================================================
double FermionFlow_1pt_Function::print_vev(int inoise, double tt,
                                           Field_F_1spinor& xi2, Field_F& xi,
                                           Field_G& U)
{
  const int    Nvol   = CommonParameters::Nvol();
  const int    NPE    = CommonParameters::NPE();
  const double ovLvol = 1.0 / Nvol / NPE;

  const int Ndim = CommonParameters::Ndim();

  const int Nquark  = m_fprop_lex.size();
  const int Nquark2 = m_params_clover->size();

  assert(Nquark == Nquark2);

  Field_F x1;
  Field_F z1, z2, z3;
  // For covariant derivative
  vector<Field_F> znu(Ndim);

  GammaMatrix *gm  = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *gm5 = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * Ndim] = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * Ndim] = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * Ndim] = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * Ndim] = gmset->get_GM(gmset->GAMMA4);
    gm5[i]           = gmset->get_GM(gmset->GAMMA5);
    delete gmset;
  }

  std::vector<dcomplex> scalar(Nquark, 0.0);
  std::vector<dcomplex> xidxi(Nquark, 0.0);
  std::vector<dcomplex> pseudo(Nquark, 0.0);
  std::vector<dcomplex> vector(Ndim * Nquark, 0.0);
  std::vector<dcomplex> axial(Ndim * Nquark, 0.0);
  std::vector<dcomplex> tensor(Ndim * Ndim * Nquark, 0.0);
  for (int id = 0; id < m_nd; id++) {
    for (int ic = 0; ic < m_nc; ++ic) {
      x1.set(0.0);
      for (int isite = 0; isite < Nvol; isite++) {
        for (int icc = 0; icc < m_nc; icc++) {
          x1.set_ri(icc, id, isite, 0, xi2.cmp_r(icc, isite, ic), xi2.cmp_i(icc, isite, ic));
        }
      }
      for (int iq = 0; iq < Nquark; ++iq) {
        int lex = iq + Nquark * (ic + m_nc * id);
        for (int isite = 0; isite < Nvol; isite++) {
          for (int idd = 0; idd < m_nd; ++idd) {
            for (int icc = 0; icc < m_nc; icc++) {
              z1.set_ri(icc, idd, isite, 0, xi.cmp_r(icc, idd, isite, lex), xi.cmp_i(icc, idd, isite, lex));
            }
          }
        }
        for (int nu = 0; nu < Ndim; ++nu) {
          m_fflow.del_symmetric(nu, znu[nu], z1, U);
        }
        dcomplex tr = dotc(x1, z1);
        scalar[iq] += tr;
        tr          = dotc(x1, x1);
        xidxi[iq]  += tr;
        mult_GM(z2, gm5[iq], z1);
        tr          = dotc(x1, z2);
        pseudo[iq] += tr;
        for (int mu = 0; mu < Ndim; ++mu) {
          int l = mu + Ndim * iq;
          mult_GM(z2, gm[mu + Ndim * iq], z1);
          tr         = dotc(x1, z2);
          vector[l] += tr;
          mult_GM(z3, gm5[iq], z2);
          tr        = dotc(x1, z3);
          axial[l] -= tr;
          for (int nu = 0; nu < Ndim; ++nu) {
            l = mu + Ndim * (nu + Ndim * iq);
            mult_GM(z3, gm[mu + Ndim * iq], znu[nu]);
            tr         = dotc(x1, z3);
            tensor[l] += tr;
          }
        }
      } // iq
    }   // ic
  }     // id

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    m_hopping_parameter = (*m_params_clover)[iq].get_double("hopping_parameter");
    scalar[iq]         *= 2.0 * m_hopping_parameter * ovLvol;
    vout.general(m_vl, " dxiSxi_sum%d %d %lf %0.16e %0.16e\n", iq, inoise, tt, real(scalar[iq]), imag(scalar[iq]));
    xidxi[iq] *= ovLvol;
    vout.general(m_vl, " dxixi_sum%d %d %lf %0.16e %0.16e\n", iq, inoise, tt, real(xidxi[iq]), imag(xidxi[iq]));
    pseudo[iq] *= 2.0 * m_hopping_parameter * ovLvol;
    vout.general(m_vl, " dxig5Sxi_sum%d %d %lf %0.16e %0.16e\n", iq, inoise, tt, real(pseudo[iq]), imag(pseudo[iq]));
    for (int mu = 0; mu < Ndim; ++mu) {
      int l = mu + Ndim * iq;
      vector[l] *= 2.0 * m_hopping_parameter * ovLvol;
      vout.general(m_vl, " dxigmuSxi_sum%d %d %d %lf %0.16e %0.16e\n", iq, inoise, mu, tt, real(vector[l]), imag(vector[l]));
      axial[l] *= 2.0 * m_hopping_parameter * ovLvol;
      vout.general(m_vl, " dxigmug5Sxi_sum%d %d %d %lf %0.16e %0.16e\n", iq, inoise, mu, tt, real(axial[l]), imag(axial[l]));
      for (int nu = 0; nu < Ndim; ++nu) {
        l          = mu + Ndim * (nu + Ndim * iq);
        tensor[l] *= 2.0 * m_hopping_parameter * ovLvol;
        vout.general(m_vl, " dpsiSxi_sum%d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, mu, nu, tt, real(tensor[l]), imag(tensor[l]));
      }
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  delete [] gm;
  delete [] gm5;
  return real(tensor[0]);
}


//====================================================================
double FermionFlow_1pt_Function::print_correlator_t(int inoise, double tt,
                                                    Field_F_1spinor& xi2,
                                                    Field_F& xi, Field_G& U)
{
  const int    Nvol          = CommonParameters::Nvol();
  const int    NPE           = CommonParameters::NPE();
  const int    Lt            = CommonParameters::Lt();
  const double normalization = double(Lt) / Nvol / NPE;

  const int Ndim = CommonParameters::Ndim();

  const int Nquark  = m_fprop_lex.size();
  const int Nquark2 = m_params_clover->size();

  assert(Nquark == Nquark2);

  Field_F x1;
  Field_F z1;
  // For covariant derivative
  vector<Field_F> znu(Ndim);

  GammaMatrix *gm   = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *g5mu = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *gm5  = new GammaMatrix[Nquark];
  GammaMatrix *gm0  = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * Ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * Ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * Ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * Ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * Ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * Ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * Ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * Ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i]             = gmset->get_GM(gmset->GAMMA5);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }

  std::vector<dcomplex> scalar(Lt * Nquark, 0.0);
  std::vector<dcomplex> xidxi(Lt * Nquark, 0.0);
  std::vector<dcomplex> pseudo(Lt * Nquark, 0.0);
  std::vector<dcomplex> vector(Lt * Ndim * Nquark, 0.0);
  std::vector<dcomplex> axial(Lt * Ndim * Nquark, 0.0);
  std::vector<dcomplex> tensor(Lt * Ndim * Ndim * Nquark, 0.0);
  for (int id = 0; id < m_nd; id++) {
    for (int ic = 0; ic < m_nc; ++ic) {
      x1.set(0.0);
      for (int isite = 0; isite < Nvol; isite++) {
        for (int icc = 0; icc < m_nc; icc++) {
          x1.set_ri(icc, id, isite, 0, xi2.cmp_r(icc, isite, ic), xi2.cmp_i(icc, isite, ic));
        }
      }
      for (int iq = 0; iq < Nquark; ++iq) {
        int lex = iq + Nquark * (ic + m_nc * id);
        for (int isite = 0; isite < Nvol; isite++) {
          for (int idd = 0; idd < m_nd; ++idd) {
            for (int icc = 0; icc < m_nc; icc++) {
              z1.set_ri(icc, idd, isite, 0, xi.cmp_r(icc, idd, isite, lex), xi.cmp_i(icc, idd, isite, lex));
            }
          }
        }
        for (int nu = 0; nu < Ndim; ++nu) {
          m_fflow.del_symmetric(nu, znu[nu], z1, U);
        }
        std::vector<dcomplex> corr(Lt);
        contract_at_t(corr, gm0[iq], z1, x1);
        for (int i = 0; i < Lt; ++i) {
          int l = i + Lt * iq;
          scalar[l] += corr[i];
        }
        contract_at_t(corr, gm0[iq], x1, x1);
        for (int i = 0; i < Lt; ++i) {
          int l = i + Lt * iq;
          xidxi[l] += corr[i];
        }
        contract_at_t(corr, gm5[iq], z1, x1);
        for (int i = 0; i < Lt; ++i) {
          int l = i + Lt * iq;
          pseudo[l] += corr[i];
        }
        for (int mu = 0; mu < Ndim; ++mu) {
          contract_at_t(corr, gm[mu + Ndim * iq], z1, x1);
          for (int i = 0; i < Lt; ++i) {
            int l = mu + Ndim * (i + Lt * iq);
            vector[l] += corr[i];
          }
          contract_at_t(corr, g5mu[mu + Ndim * iq], z1, x1);
          for (int i = 0; i < Lt; ++i) {
            int l = mu + Ndim * (i + Lt * iq);
            axial[l] -= corr[i];
          }
          for (int nu = 0; nu < Ndim; ++nu) {
            contract_at_t(corr, gm[mu + Ndim * iq], znu[nu], x1);
            for (int i = 0; i < Lt; ++i) {
              int l = mu + Ndim * (nu + Ndim * (i + Lt * iq));
              tensor[l] += corr[i];
            }
          }
        }
      } // iq
    }   // ic
  }     // id

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    m_hopping_parameter = (*m_params_clover)[iq].get_double("hopping_parameter");
    dcomplex dxiSxi   = 0;
    dcomplex dxixi    = 0;
    dcomplex dxig5Sxi = 0;
    for (int i = 0; i < Lt; ++i) {
      int      lex = i + Lt * iq;
      dcomplex scr = scalar[lex] * 2.0 * m_hopping_parameter * normalization;
      vout.general(m_vl, " dxiSxi_t%d= %d %lf %d %0.16e %0.16e\n", iq, inoise, tt, i, real(scr), imag(scr));
      dxiSxi += scr;
      scr     = xidxi[lex] * normalization;
      vout.general(m_vl, " dxixi_t%d= %d %lf %d %0.16e %0.16e\n", iq, inoise, tt, i, real(scr), imag(scr));
      dxixi += scr;
      scr    = pseudo[lex] * 2.0 * m_hopping_parameter * normalization;
      vout.general(m_vl, " dxig5Sxi_t%d= %d %lf %d %0.16e %0.16e\n", iq, inoise, tt, i, real(scr), imag(scr));
      dxig5Sxi += scr;
    }
    vout.general(m_vl, " dxiSxi_sumt%d= %d %lf %0.16e %0.16e\n", iq, inoise, tt, real(dxiSxi) / Lt, imag(dxiSxi) / Lt);
    vout.general(m_vl, " dxixi_sumt%d= %d %lf %0.16e %0.16e\n", iq, inoise, tt, real(dxixi) / Lt, imag(dxixi) / Lt);
    vout.general(m_vl, " dxig5Sxi_sumt%d= %d %lf %0.16e %0.16e\n", iq, inoise, tt, real(dxig5Sxi) / Lt, imag(dxig5Sxi) / Lt);
    for (int mu = 0; mu < Ndim; ++mu) {
      dcomplex dxigmuSxi   = 0;
      dcomplex dxigmug5Sxi = 0;
      for (int i = 0; i < Lt; ++i) {
        int      lex = mu + Ndim * (i + Lt * iq);
        dcomplex scr = vector[lex] * 2.0 * m_hopping_parameter * normalization;
        vout.general(m_vl, " dxigmuSxi_t%d= %d %d %lf %d %0.16e %0.16e\n", iq, inoise, mu, tt, i, real(scr), imag(scr));
        dxigmuSxi += scr;
        scr        = axial[lex] * 2.0 * m_hopping_parameter * normalization;
        vout.general(m_vl, " dxigmug5Sxi_t%d= %d %d %lf %d %0.16e %0.16e\n", iq, inoise, mu, tt, i, real(scr), imag(scr));
        dxigmug5Sxi += scr;
      }
      vout.general(m_vl, " dxigmuSxi_sumt%d= %d %d %lf %0.16e %0.16e\n", iq, inoise, mu, tt, real(dxigmuSxi) / Lt, imag(dxigmuSxi) / Lt);
      vout.general(m_vl, " dxigmug5Sxi_sumt%d= %d %d %lf %0.16e %0.16e\n", iq, inoise, mu, tt, real(dxigmug5Sxi) / Lt, imag(dxigmug5Sxi) / Lt);
    }
    for (int mu = 0; mu < Ndim; ++mu) {
      for (int nu = 0; nu < Ndim; ++nu) {
        dcomplex dpsiSxi = 0;
        for (int i = 0; i < Lt; ++i) {
          int      lex = mu + Ndim * (nu + Ndim * (i + Lt * iq));
          dcomplex scr = tensor[lex] * 2.0 * m_hopping_parameter * normalization;
          vout.general(m_vl, " dpsiSxi_t%d= %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, mu, nu, tt, i, real(scr), imag(scr));
          dpsiSxi += scr;
        }
        vout.general(m_vl, " dpsiSxi_sumt%d= %d %d %d %lf %0.16e %0.16e\n", iq, inoise, mu, nu, tt, real(dpsiSxi) / Lt, imag(dpsiSxi) / Lt);
      }
    }
  }
  double result = real(tensor[0] * 2.0 * m_hopping_parameter * normalization);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  return result;
}


//====================================================================
double FermionFlow_1pt_Function::print_correlator_t_FT(int inoise, double tt,
                                                       Field_F_1spinor& xi2,
                                                       Field_F& xi, Field_G& U)
{
  const int    Nvol          = CommonParameters::Nvol();
  const int    NPE           = CommonParameters::NPE();
  const int    Lt            = CommonParameters::Lt();
  const double normalization = double(Lt) / Nvol / NPE;

  const int Ndim = CommonParameters::Ndim();

  double result = 0.0;

  const int Np = (2 * m_max_mom + 1);

  vector<int> source_position(4, 0);
  vector<int> momentum_sink(3);

  const int Nquark  = m_fprop_lex.size();
  const int Nquark2 = m_params_clover->size();
  assert(Nquark == Nquark2);

  Field_F x1;
  Field_F z1;
  // For covariant derivative
  vector<Field_F> znu(Ndim);

  GammaMatrix *gm   = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *g5mu = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *gm5  = new GammaMatrix[Nquark];
  GammaMatrix *gm0  = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * Ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * Ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * Ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * Ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * Ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * Ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * Ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * Ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i]             = gmset->get_GM(gmset->GAMMA5);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }

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
        std::vector<dcomplex> scalar(Lt * Nquark, 0.0);
        std::vector<dcomplex> xidxi(Lt * Nquark, 0.0);
        std::vector<dcomplex> pseudo(Lt * Nquark, 0.0);
        std::vector<dcomplex> vector(Lt * Ndim * Nquark, 0.0);
        std::vector<dcomplex> axial(Lt * Ndim * Nquark, 0.0);
        std::vector<dcomplex> tensor(Lt * Ndim * Ndim * Nquark, 0.0);
        for (int id = 0; id < m_nd; id++) {
          for (int ic = 0; ic < m_nc; ++ic) {
            x1.set(0.0);
            for (int isite = 0; isite < Nvol; isite++) {
              for (int icc = 0; icc < m_nc; icc++) {
                x1.set_ri(icc, id, isite, 0, xi2.cmp_r(icc, isite, ic), xi2.cmp_i(icc, isite, ic));
              }
            }
            for (int iq = 0; iq < Nquark; ++iq) {
              int lex = iq + Nquark * (ic + m_nc * id);
              for (int isite = 0; isite < Nvol; isite++) {
                for (int idd = 0; idd < m_nd; ++idd) {
                  for (int icc = 0; icc < m_nc; icc++) {
                    z1.set_ri(icc, idd, isite, 0, xi.cmp_r(icc, idd, isite, lex), xi.cmp_i(icc, idd, isite, lex));
                  }
                }
              }
              for (int nu = 0; nu < Ndim; ++nu) {
                m_fflow.del_symmetric(nu, znu[nu], z1, U);
              }
              std::vector<dcomplex> corr(Lt);
              contract_at_t(corr, momentum_sink, gm0[iq], source_position, z1, x1);
              for (int i = 0; i < Lt; ++i) {
                int l = i + Lt * iq;
                scalar[l] += corr[i];
              }
              contract_at_t(corr, momentum_sink, gm0[iq], source_position, x1, x1);
              for (int i = 0; i < Lt; ++i) {
                int l = i + Lt * iq;
                xidxi[l] += corr[i];
              }
              contract_at_t(corr, momentum_sink, gm5[iq], source_position, z1, x1);
              for (int i = 0; i < Lt; ++i) {
                int l = i + Lt * iq;
                pseudo[l] += corr[i];
              }
              for (int mu = 0; mu < Ndim; ++mu) {
                contract_at_t(corr, momentum_sink, gm[mu + Ndim * iq], source_position, z1, x1);
                for (int i = 0; i < Lt; ++i) {
                  int l = mu + Ndim * (i + Lt * iq);
                  vector[l] += corr[i];
                }
                contract_at_t(corr, momentum_sink, g5mu[mu + Ndim * iq], source_position, z1, x1);
                for (int i = 0; i < Lt; ++i) {
                  int l = mu + Ndim * (i + Lt * iq);
                  axial[l] -= corr[i];
                }
                for (int nu = 0; nu < Ndim; ++nu) {
                  contract_at_t_cos(corr, momentum_sink, gm[mu + Ndim * iq], source_position, znu[nu], x1);
                  for (int i = 0; i < Lt; ++i) {
                    int l = mu + Ndim * (nu + Ndim * (i + Lt * iq));
                    tensor[l] += corr[i];
                  }
                }
              }
            } //iq
          }   // ic
        }     // id

        for (int iq = 0; iq < Nquark; ++iq) {
          m_hopping_parameter = (*m_params_clover)[iq].get_double("hopping_parameter");
          dcomplex dxiSxi   = 0;
          dcomplex dxixi    = 0;
          dcomplex dxig5Sxi = 0;
          for (int i = 0; i < Lt; ++i) {
            int      lex = i + Lt * iq;
            dcomplex scr = scalar[lex] * 2.0 * m_hopping_parameter * normalization;
            vout.general(m_vl, " dxiSxi_t_FT%d= %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, i, real(scr), imag(scr));
            dxiSxi += scr;
            scr     = xidxi[lex] * normalization;
            vout.general(m_vl, " dxixi_t_FT%d= %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, i, real(scr), imag(scr));
            dxixi += scr;
            scr    = pseudo[lex] * 2.0 * m_hopping_parameter * normalization;
            vout.general(m_vl, " dxig5Sxi_t_FT%d= %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, i, real(scr), imag(scr));
            dxig5Sxi += scr;
          }
          vout.general(m_vl, " dxiSxi_sumt_FT%d= %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, real(dxiSxi) / Lt, imag(dxiSxi) / Lt);
          vout.general(m_vl, " dxixi_sumt_FT%d= %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, real(dxixi) / Lt, imag(dxixi) / Lt);
          vout.general(m_vl, " dxig5Sxi_sumt_FT%d= %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, real(dxig5Sxi) / Lt, imag(dxig5Sxi) / Lt);
          for (int mu = 0; mu < Ndim; ++mu) {
            dcomplex dxigmuSxi   = 0;
            dcomplex dxigmug5Sxi = 0;
            for (int i = 0; i < Lt; ++i) {
              int      lex = mu + Ndim * (i + Lt * iq);
              dcomplex scr = vector[lex] * 2.0 * m_hopping_parameter * normalization;
              vout.general(m_vl, " dxigmuSxi_t_FT%d= %d %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, i, real(scr), imag(scr));
              dxigmuSxi += scr;
              scr        = axial[lex] * 2.0 * m_hopping_parameter * normalization;
              vout.general(m_vl, " dxigmug5Sxi_t_FT%d= %d %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, i, real(scr), imag(scr));
              dxigmug5Sxi += scr;
            }
            vout.general(m_vl, " dxigmuSxi_sumt_FT%d= %d %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, real(dxigmuSxi) / Lt, imag(dxigmuSxi) / Lt);
            vout.general(m_vl, " dxigmug5Sxi_sumt_FT%d= %d %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, real(dxigmug5Sxi) / Lt, imag(dxigmug5Sxi) / Lt);
          }
          for (int mu = 0; mu < Ndim; ++mu) {
            for (int nu = 0; nu < Ndim; ++nu) {
              dcomplex dpsiSxi = 0;
              for (int i = 0; i < Lt; ++i) {
                int      lex = mu + Ndim * (nu + Ndim * (i + Lt * iq));
                dcomplex scr = tensor[lex] * 2.0 * m_hopping_parameter * normalization;
                vout.general(m_vl, " dpsiSxi_t_FT%d= %d %d %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, i, real(scr), imag(scr));
                result   = real(scr);
                dpsiSxi += scr;
              }
              vout.general(m_vl, " dpsiSxi_sumt_FT%d= %d %d %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, real(dpsiSxi) / Lt, imag(dpsiSxi) / Lt);
            }
          }
        } //iq
      }   // ipz
    }     // ipy
  }       // ipx
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  return result;
}


//====================================================================
double FermionFlow_1pt_Function::print_correlator_x(int inoise, double tt,
                                                    Field_F_1spinor& xi2,
                                                    Field_F& xi, Field_G& U)
{
  const int    Nvol          = CommonParameters::Nvol();
  const int    NPE           = CommonParameters::NPE();
  const int    Lx            = CommonParameters::Lx();
  const double normalization = double(Lx) / Nvol / NPE;

  const int Ndim = CommonParameters::Ndim();

  const int Nquark  = m_fprop_lex.size();
  const int Nquark2 = m_params_clover->size();

  assert(Nquark == Nquark2);

  Field_F x1;
  Field_F z1;
  // For covariant derivative
  vector<Field_F> znu(Ndim);

  GammaMatrix *gm   = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *g5mu = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *gm5  = new GammaMatrix[Nquark];
  GammaMatrix *gm0  = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * Ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * Ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * Ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * Ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * Ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * Ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * Ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * Ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i]             = gmset->get_GM(gmset->GAMMA5);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }

  std::vector<dcomplex> scalar(Lx * Nquark, 0.0);
  std::vector<dcomplex> xidxi(Lx * Nquark, 0.0);
  std::vector<dcomplex> pseudo(Lx * Nquark, 0.0);
  std::vector<dcomplex> vector(Lx * Ndim * Nquark, 0.0);
  std::vector<dcomplex> axial(Lx * Ndim * Nquark, 0.0);
  std::vector<dcomplex> tensor(Lx * Ndim * Ndim * Nquark, 0.0);
  for (int id = 0; id < m_nd; id++) {
    for (int ic = 0; ic < m_nc; ++ic) {
      x1.set(0.0);
      for (int isite = 0; isite < Nvol; isite++) {
        for (int icc = 0; icc < m_nc; icc++) {
          x1.set_ri(icc, id, isite, 0, xi2.cmp_r(icc, isite, ic), xi2.cmp_i(icc, isite, ic));
        }
      }
      for (int iq = 0; iq < Nquark; ++iq) {
        int lex = iq + Nquark * (ic + m_nc * id);
        for (int isite = 0; isite < Nvol; isite++) {
          for (int idd = 0; idd < m_nd; ++idd) {
            for (int icc = 0; icc < m_nc; icc++) {
              z1.set_ri(icc, idd, isite, 0, xi.cmp_r(icc, idd, isite, lex), xi.cmp_i(icc, idd, isite, lex));
            }
          }
        }
        for (int nu = 0; nu < Ndim; ++nu) {
          m_fflow.del_symmetric(nu, znu[nu], z1, U);
        }
        std::vector<dcomplex> corr(Lx);
        contract_at_x(corr, gm0[iq], z1, x1);
        for (int i = 0; i < Lx; ++i) {
          int l = i + Lx * iq;
          scalar[l] += corr[i];
        }
        contract_at_x(corr, gm0[iq], x1, x1);
        for (int i = 0; i < Lx; ++i) {
          int l = i + Lx * iq;
          xidxi[l] += corr[i];
        }
        contract_at_x(corr, gm5[iq], z1, x1);
        for (int i = 0; i < Lx; ++i) {
          int l = i + Lx * iq;
          pseudo[l] += corr[i];
        }
        for (int mu = 0; mu < Ndim; ++mu) {
          contract_at_x(corr, gm[mu + Ndim * iq], z1, x1);
          for (int i = 0; i < Lx; ++i) {
            int l = mu + Ndim * (i + Lx * iq);
            vector[l] += corr[i];
          }
          contract_at_x(corr, g5mu[mu + Ndim * iq], z1, x1);
          for (int i = 0; i < Lx; ++i) {
            int l = mu + Ndim * (i + Lx * iq);
            axial[l] -= corr[i];
          }
          for (int nu = 0; nu < Ndim; ++nu) {
            contract_at_x(corr, gm[mu + Ndim * iq], znu[nu], x1);
            for (int i = 0; i < Lx; ++i) {
              int l = mu + Ndim * (nu + Ndim * (i + Lx * iq));
              tensor[l] += corr[i];
            }
          }
        }
      } // iq
    }   // ic
  }     // id

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int iq = 0; iq < Nquark; ++iq) {
    m_hopping_parameter = (*m_params_clover)[iq].get_double("hopping_parameter");
    dcomplex dxiSxi   = 0;
    dcomplex dxixi    = 0;
    dcomplex dxig5Sxi = 0;
    for (int i = 0; i < Lx; ++i) {
      int      lex = i + Lx * iq;
      dcomplex scr = scalar[lex] * 2.0 * m_hopping_parameter * normalization;
      vout.general(m_vl, " dxiSxi_x%d= %d %lf %d %0.16e %0.16e\n", iq, inoise, tt, i, real(scr), imag(scr));
      dxiSxi += scr;
      scr     = xidxi[lex] * normalization;
      vout.general(m_vl, " dxixi_x%d= %d %lf %d %0.16e %0.16e\n", iq, inoise, tt, i, real(scr), imag(scr));
      dxixi += scr;
      scr    = pseudo[lex] * 2.0 * m_hopping_parameter * normalization;
      vout.general(m_vl, " dxig5Sxi_x%d= %d %lf %d %0.16e %0.16e\n", iq, inoise, tt, i, real(scr), imag(scr));
      dxig5Sxi += scr;
    }
    vout.general(m_vl, " dxiSxi_sumx%d= %d %lf %0.16e %0.16e\n", iq, inoise, tt, real(dxiSxi) / Lx, imag(dxiSxi) / Lx);
    vout.general(m_vl, " dxixi_sumx%d= %d %lf %0.16e %0.16e\n", iq, inoise, tt, real(dxixi) / Lx, imag(dxixi) / Lx);
    vout.general(m_vl, " dxig5Sxi_sumx%d= %d %lf %0.16e %0.16e\n", iq, inoise, tt, real(dxig5Sxi) / Lx, imag(dxig5Sxi) / Lx);
    for (int mu = 0; mu < Ndim; ++mu) {
      dcomplex dxigmuSxi   = 0;
      dcomplex dxigmug5Sxi = 0;
      for (int i = 0; i < Lx; ++i) {
        int      lex = mu + Ndim * (i + Lx * iq);
        dcomplex scr = vector[lex] * 2.0 * m_hopping_parameter * normalization;
        vout.general(m_vl, " dxigmuSxi_x%d= %d %d %lf %d %0.16e %0.16e\n", iq, inoise, mu, tt, i, real(scr), imag(scr));
        dxigmuSxi += scr;
        scr        = axial[lex] * 2.0 * m_hopping_parameter * normalization;
        vout.general(m_vl, " dxigmug5Sxi_x%d= %d %d %lf %d %0.16e %0.16e\n", iq, inoise, mu, tt, i, real(scr), imag(scr));
        dxigmug5Sxi += scr;
      }
      vout.general(m_vl, " dxigmuSxi_sumx%d= %d %d %lf %0.16e %0.16e\n", iq, inoise, mu, tt, real(dxigmuSxi) / Lx, imag(dxigmuSxi) / Lx);
      vout.general(m_vl, " dxigmug5Sxi_sumx%d= %d %d %lf %0.16e %0.16e\n", iq, inoise, mu, tt, real(dxigmug5Sxi) / Lx, imag(dxigmug5Sxi) / Lx);
    }
    for (int mu = 0; mu < Ndim; ++mu) {
      for (int nu = 0; nu < Ndim; ++nu) {
        dcomplex dpsiSxi = 0;
        for (int i = 0; i < Lx; ++i) {
          int      lex = mu + Ndim * (nu + Ndim * (i + Lx * iq));
          dcomplex scr = tensor[lex] * 2.0 * m_hopping_parameter * normalization;
          vout.general(m_vl, " dpsiSxi_x%d= %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, mu, nu, tt, i, real(scr), imag(scr));
          dpsiSxi += scr;
        }
        vout.general(m_vl, " dpsiSxi_sumx%d= %d %d %d %lf %0.16e %0.16e\n", iq, inoise, mu, nu, tt, real(dpsiSxi) / Lx, imag(dpsiSxi) / Lx);
      }
    }
  }
  double result = real(tensor[0] * 2.0 * m_hopping_parameter * normalization);
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  return result;
}


//====================================================================
double FermionFlow_1pt_Function::print_correlator_x_FT(int inoise, double tt,
                                                       Field_F_1spinor& xi2,
                                                       Field_F& xi, Field_G& U)
{
  const int    Nvol          = CommonParameters::Nvol();
  const int    NPE           = CommonParameters::NPE();
  const int    Lx            = CommonParameters::Lx();
  const double normalization = double(Lx) / Nvol / NPE;

  const int Ndim = CommonParameters::Ndim();

  double result = 0.0;

  const int Np = (2 * m_max_mom + 1);

  vector<int> source_position(Ndim, 0);
  vector<int> momentum_sink(Ndim - 1);

  const int Nquark  = m_fprop_lex.size();
  const int Nquark2 = m_params_clover->size();
  assert(Nquark == Nquark2);

  Field_F x1;
  Field_F z1;
  // For covariant derivative
  vector<Field_F> znu(Ndim);

  GammaMatrix *gm   = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *g5mu = new GammaMatrix[Ndim * Nquark];
  GammaMatrix *gm5  = new GammaMatrix[Nquark];
  GammaMatrix *gm0  = new GammaMatrix[Nquark];
  for (int i = 0; i < Nquark; ++i) {
    m_str_gmset_type = (*m_params_clover)[i].get_string("gamma_matrix_type");
    GammaMatrixSet *gmset = GammaMatrixSet::New(m_str_gmset_type);
    gm[0 + i * Ndim]   = gmset->get_GM(gmset->GAMMA1);
    gm[1 + i * Ndim]   = gmset->get_GM(gmset->GAMMA2);
    gm[2 + i * Ndim]   = gmset->get_GM(gmset->GAMMA3);
    gm[3 + i * Ndim]   = gmset->get_GM(gmset->GAMMA4);
    g5mu[0 + i * Ndim] = gmset->get_GM(gmset->GAMMA51);
    g5mu[1 + i * Ndim] = gmset->get_GM(gmset->GAMMA52);
    g5mu[2 + i * Ndim] = gmset->get_GM(gmset->GAMMA53);
    g5mu[3 + i * Ndim] = gmset->get_GM(gmset->GAMMA54);
    gm5[i]             = gmset->get_GM(gmset->GAMMA5);
    gm0[i]             = gmset->get_GM(gmset->UNITY);
    delete gmset;
  }

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
        std::vector<dcomplex> scalar(Lx * Nquark, 0.0);
        std::vector<dcomplex> xidxi(Lx * Nquark, 0.0);
        std::vector<dcomplex> pseudo(Lx * Nquark, 0.0);
        std::vector<dcomplex> vector(Lx * Ndim * Nquark, 0.0);
        std::vector<dcomplex> axial(Lx * Ndim * Nquark, 0.0);
        std::vector<dcomplex> tensor(Lx * Ndim * Ndim * Nquark, 0.0);
        for (int id = 0; id < m_nd; id++) {
          for (int ic = 0; ic < m_nc; ++ic) {
            x1.set(0.0);
            for (int isite = 0; isite < Nvol; isite++) {
              for (int icc = 0; icc < m_nc; icc++) {
                x1.set_ri(icc, id, isite, 0, xi2.cmp_r(icc, isite, ic), xi2.cmp_i(icc, isite, ic));
              }
            }
            for (int iq = 0; iq < Nquark; ++iq) {
              int lex = iq + Nquark * (ic + m_nc * id);
              for (int isite = 0; isite < Nvol; isite++) {
                for (int idd = 0; idd < m_nd; ++idd) {
                  for (int icc = 0; icc < m_nc; icc++) {
                    z1.set_ri(icc, idd, isite, 0, xi.cmp_r(icc, idd, isite, lex), xi.cmp_i(icc, idd, isite, lex));
                  }
                }
              }
              for (int nu = 0; nu < Ndim; ++nu) {
                m_fflow.del_symmetric(nu, znu[nu], z1, U);
              }
              std::vector<dcomplex> corr(Lx);
              contract_at_x(corr, momentum_sink, gm0[iq], source_position, z1, x1);
              for (int i = 0; i < Lx; ++i) {
                int l = i + Lx * iq;
                scalar[l] += corr[i];
              }
              contract_at_x(corr, momentum_sink, gm0[iq], source_position, x1, x1);
              for (int i = 0; i < Lx; ++i) {
                int l = i + Lx * iq;
                xidxi[l] += corr[i];
              }
              contract_at_x(corr, momentum_sink, gm5[iq], source_position, z1, x1);
              for (int i = 0; i < Lx; ++i) {
                int l = i + Lx * iq;
                pseudo[l] += corr[i];
              }
              for (int mu = 0; mu < Ndim; ++mu) {
                contract_at_x(corr, momentum_sink, gm[mu + Ndim * iq], source_position, z1, x1);
                for (int i = 0; i < Lx; ++i) {
                  int l = mu + Ndim * (i + Lx * iq);
                  vector[l] += corr[i];
                }
                contract_at_x(corr, momentum_sink, g5mu[mu + Ndim * iq], source_position, z1, x1);
                for (int i = 0; i < Lx; ++i) {
                  int l = mu + Ndim * (i + Lx * iq);
                  axial[l] -= corr[i];
                }
                for (int nu = 0; nu < Ndim; ++nu) {
                  contract_at_x_cos(corr, momentum_sink, gm[mu + Ndim * iq], source_position, znu[nu], x1);
                  for (int i = 0; i < Lx; ++i) {
                    int l = mu + Ndim * (nu + Ndim * (i + Lx * iq));
                    tensor[l] += corr[i];
                  }
                }
              }
            } //iq
          }   // ic
        }     // id

        for (int iq = 0; iq < Nquark; ++iq) {
          m_hopping_parameter = (*m_params_clover)[iq].get_double("hopping_parameter");
          dcomplex dxiSxi   = 0;
          dcomplex dxixi    = 0;
          dcomplex dxig5Sxi = 0;
          for (int i = 0; i < Lx; ++i) {
            int      lex = i + Lx * iq;
            dcomplex scr = scalar[lex] * 2.0 * m_hopping_parameter * normalization;
            vout.general(m_vl, " dxiSxi_x_FT%d= %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, i, real(scr), imag(scr));
            dxiSxi += scr;
            scr     = xidxi[lex] * normalization;
            vout.general(m_vl, " dxixi_x_FT%d= %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, i, real(scr), imag(scr));
            dxixi += scr;
            scr    = pseudo[lex] * 2.0 * m_hopping_parameter * normalization;
            vout.general(m_vl, " dxig5Sxi_x_FT%d= %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, i, real(scr), imag(scr));
            dxig5Sxi += scr;
          }
          vout.general(m_vl, " dxiSxi_sumx_FT%d= %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, real(dxiSxi) / Lx, imag(dxiSxi) / Lx);
          vout.general(m_vl, " dxixi_sumx_FT%d= %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, real(dxixi) / Lx, imag(dxixi) / Lx);
          vout.general(m_vl, " dxig5Sxi_sumx_FT%d= %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, real(dxig5Sxi) / Lx, imag(dxig5Sxi) / Lx);
          for (int mu = 0; mu < Ndim; ++mu) {
            dcomplex dxigmuSxi   = 0;
            dcomplex dxigmug5Sxi = 0;
            for (int i = 0; i < Lx; ++i) {
              int      lex = mu + Ndim * (i + Lx * iq);
              dcomplex scr = vector[lex] * 2.0 * m_hopping_parameter * normalization;
              vout.general(m_vl, " dxigmuSxi_x_FT%d= %d %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, i, real(scr), imag(scr));
              dxigmuSxi += scr;
              scr        = axial[lex] * 2.0 * m_hopping_parameter * normalization;
              vout.general(m_vl, " dxigmug5Sxi_x_FT%d= %d %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, i, real(scr), imag(scr));
              dxigmug5Sxi += scr;
            }
            vout.general(m_vl, " dxigmuSxi_sumx_FT%d= %d %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, real(dxigmuSxi) / Lx, imag(dxigmuSxi) / Lx);
            vout.general(m_vl, " dxigmug5Sxi_sumx_FT%d= %d %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, tt, real(dxigmug5Sxi) / Lx, imag(dxigmug5Sxi) / Lx);
          }
          for (int mu = 0; mu < Ndim; ++mu) {
            for (int nu = 0; nu < Ndim; ++nu) {
              dcomplex dpsiSxi = 0;
              for (int i = 0; i < Lx; ++i) {
                int      lex = mu + Ndim * (nu + Ndim * (i + Lx * iq));
                dcomplex scr = tensor[lex] * 2.0 * m_hopping_parameter * normalization;
                vout.general(m_vl, " dpsiSxi_x_FT%d= %d %d %d %d %d %d %lf %d %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, i, real(scr), imag(scr));
                result   = real(scr);
                dpsiSxi += scr;
              }
              vout.general(m_vl, " dpsiSxi_sumx_FT%d= %d %d %d %d %d %d %lf %0.16e %0.16e\n", iq, inoise, momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, tt, real(dpsiSxi) / Lx, imag(dpsiSxi) / Lx);
            }
          }
        } //iq
      }   // ipz
    }     // ipy
  }       // ipx
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  delete [] gm;
  delete [] g5mu;
  delete [] gm5;
  delete [] gm0;
  return result;
}


//====================================================================
//============================================================END=====
