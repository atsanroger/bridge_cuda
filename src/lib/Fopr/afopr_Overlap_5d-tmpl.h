/*!
        @file    afopr_Overlap_5d-tmpl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-05 20:48:53 #$

        @version $LastChangedRevision: 2530 $
*/

#include "Fopr/afopr_Overlap_5d.h"
#include "lib/ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = AFopr_Overlap_5d<AFIELD>::register_factory();
}
#endif

template<typename AFIELD>
const std::string AFopr_Overlap_5d<AFIELD>::class_name
                                       = "AFopr_Overlap_5d<AFIELD>";

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.general(m_vl, "  -- For 5d overlap solver\n");
  vout.increase_indent();

  std::string kernel_type;
  int err = params.fetch_string("kernel_type", kernel_type);
  if(err > 0){
    vout.crucial(m_vl, "Error at %s: kernel_type is not specified.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  m_kernel_type = kernel_type;
  std::string kernel_eo_type = kernel_type + "_eo";

  double M0;
  params.fetch_double("domain_wall_height", M0);

  Parameters params_kernel = params;
  double kappa = 1.0/(8.0 - 2.0 * M0);
  params_kernel.set_double("hopping_parameter", kappa);

  m_fopr_w = AFOPR::New(kernel_eo_type, params_kernel);
  m_kernel_created  = true;

  m_Nin   = m_fopr_w->field_nin();
  m_Nvol2 = m_fopr_w->field_nvol();
  m_Nex   = 0;
  m_Ndim = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  m_Nsbt = 0;
  m_mode = "DdagD_eo";

  m_z1.reset(m_Nin, m_Nvol2, 1);
  m_z2.reset(m_Nin, m_Nvol2, 1);

  m_w1.reset(m_Nin, m_Nvol2, 1);
  m_v1.reset(m_Nin, m_Nvol2, 1);
  m_v2.reset(m_Nin, m_Nvol2, 1);

  set_parameters(params);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::init(AFOPR *fopr, const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.general(m_vl, "  -- For 5d overlap solver\n");
  vout.increase_indent();

  m_fopr_w = fopr;
  m_kernel_created = false;

  m_kernel_type = "unknown";
  m_repr        = "unknown";

  m_Nin   = m_fopr_w->field_nin();
  m_Nvol2 = m_fopr_w->field_nvol();
  m_Nex   = 0;
  m_Ndim = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  m_Nsbt = 0;
  m_mode = "DdagD_eo";

  m_z1.reset(m_Nin, m_Nvol2, 1);
  m_z2.reset(m_Nin, m_Nvol2, 1);

  m_w1.reset(m_Nin, m_Nvol2, 1);
  m_v1.reset(m_Nin, m_Nvol2, 1);
  m_v2.reset(m_Nin, m_Nvol2, 1);

  set_parameters(params);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n", class_name.c_str());
  vout.general(m_vl, "  -- For 5d overlap solver\n");
  vout.increase_indent();

  m_kernel_created = false;

  m_kernel_type = "unknown";
  m_repr        = "unknown";

  m_Nin   = m_fopr_w->field_nin();
  m_Nvol2 = m_fopr_w->field_nvol();
  m_Nex   = 0;
  m_Ndim = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  m_Nsbt = 0;
  m_mode = "DdagD_eo";

  m_z1.reset(m_Nin, m_Nvol2, 1);
  m_z2.reset(m_Nin, m_Nvol2, 1);

  m_w1.reset(m_Nin, m_Nvol2, 1);
  m_v1.reset(m_Nin, m_Nvol2, 1);
  m_v2.reset(m_Nin, m_Nvol2, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::tidyup()
{
  if(m_kernel_created == true) delete m_fopr_w;
}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double mq, M0;
  int    Np;
  double x_min, x_max;
  int    Niter_ms;
  double Stop_cond_ms;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);
  err += params.fetch_int("maximum_number_of_iteration", Niter_ms);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond_ms);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
	         class_name.c_str());
    exit(EXIT_FAILURE);
  }

  std::string repr, kernel;
  if (!params.fetch_string("gamma_matrix_type", repr)) m_repr = repr;
  if (!params.fetch_string("kernel_type", kernel)) m_kernel_type = kernel;

  set_parameters(real_t(mq), real_t(M0),
                 Np, real_t(x_min), real_t(x_max),
                 Niter_ms, real_t(Stop_cond_ms), bc);
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::set_parameters(const real_t mq,
                                              const real_t M0,
                                              const int Np,
                                              const real_t x_min,
                                              const real_t x_max,
                                              const int Niter_ms,
                                              const real_t Stop_cond_ms,
                                              const std::vector<int> bc)
{
#pragma omp barrier

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(mq);
  err += ParameterCheck::non_zero(M0);
  err += ParameterCheck::non_zero(Np);
  // NB. x_min,x_max == 0 is allowed.
  err += ParameterCheck::non_zero(Niter_ms);
  err += ParameterCheck::square_non_zero(Stop_cond_ms);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();

  if(ith == 0){
    m_mq = mq;
    m_M0 = M0;

    m_Np    = Np;
    m_x_min = x_min;
    m_x_max = x_max;

    m_Niter_ms     = Niter_ms;
    m_Stop_cond_ms = Stop_cond_ms;

    m_boundary = bc;
  }

  // print input parameters
  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq    = %12.8f\n", m_mq);
  vout.general(m_vl, "  M0    = %12.8f\n", m_M0);
  vout.general(m_vl, "  Np    = %4d\n", m_Np);
  vout.general(m_vl, "  x_min = %12.8f\n", m_x_min);
  vout.general(m_vl, "  x_max = %12.8f\n", m_x_max);
  vout.general(m_vl, "  Niter_ms     = %8d\n", m_Niter_ms);
  vout.general(m_vl, "  Stop_cond_ms = %8.2e\n", m_Stop_cond_ms);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

  set_coefficients();

  if(ith == 0){
    m_Nex = 2 * m_Np + 1;
    if(m_t1.nex() != m_Nex){
      m_t1.reset(m_Nin, m_Nvol2, m_Nex);
      m_t2.reset(m_Nin, m_Nvol2, m_Nex);
      m_t3.reset(m_Nin, m_Nvol2, m_Nex);
    }
  }
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::set_coefficients()
{
  int ith = ThreadManager::get_thread_id();

  if(ith == 0){

    // Zolotarev coefficients and shift values
    m_sigma.resize(m_Np);
    m_cl.resize(2 * m_Np);
    m_bl.resize(m_Np);

    // Zolotarev coefficient defined
    const real_t bmax = m_x_max / m_x_min;

    Math_Sign_Zolotarev sign_func(m_Np, bmax);
    sign_func.get_sign_parameters(m_cl, m_bl);

    for (int i = 0; i < m_Np; i++) {
      m_sigma[i] = m_cl[2 * i] * m_x_min * m_x_min;
    }

    for (int i = 0; i < m_Np; i++) {
      vout.general(m_vl, " %3d %12.4e %12.4e %12.4e\n",
                 i, m_cl[i], m_cl[i + m_Np], m_bl[i]);
    }

    m_p_sqrt.resize(m_Np);
    m_q_sqrt.resize(m_Np);

    real_t b_sum = 0.0;
    for (int ip = 0; ip < m_Np; ++ip) {
      int    ik    = m_Np - ip - 1;
      real_t p_tmp = (m_M0 - 0.5 * m_mq) * m_bl[ip] *
                     (m_cl[2 * m_Np - 1] - m_cl[2 * ip]) * m_x_min;
      m_p_sqrt[ik] = sqrt(p_tmp);
      m_q_sqrt[ik] = sqrt(m_cl[2 * ip] * m_x_min * m_x_min);
      b_sum        = b_sum + m_bl[ip];
    }

    m_p0_parameter = (m_M0 - 0.5 * m_mq) * b_sum / m_x_min;
    m_R_parameter  = m_M0 + 0.5 * m_mq;
    m_h            = 1.0;

    m_rl.resize(m_Np);
    m_sl.resize(m_Np);
    for (int ik = 0; ik < m_Np; ++ik) {
      m_rl[ik] = -m_q_sqrt[ik];
      m_sl[ik] = -m_p_sqrt[ik] / (1.0 + m_q_sqrt[ik] * m_q_sqrt[ik]);
    }

    m_u0 = 0.0;
    for (int ik = 0; ik < m_Np; ++ik) {
      m_u0 = m_u0 - m_sl[ik] * m_p_sqrt[ik];
    }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", double(m_mq));
  params.set_double("domain_wall_height", double(m_M0));
  params.set_int("number_of_poles", m_Np);
  params.set_double("lower_bound", double(m_x_min));
  params.set_double("upper_bound", double(m_x_max));
  params.set_int("maximum_number_of_iteration", m_Niter_ms);
  params.set_double("convergence_criterion_squared", double(m_Stop_cond_ms));
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("kernel_type", m_kernel_type);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::set_lowmodes(const int Nsbt,
                                            std::vector<real_t> *ev,
                                            std::vector<AFIELD> *vk)
{
  m_Nsbt = Nsbt;

  vout.general(m_vl, "  Nsbt = %d", Nsbt);

  if (m_Nsbt > 0) {
    m_ev.resize(m_Nsbt);
    m_vk.resize(2 * m_Nsbt);

    const int Nin   = (*vk)[0].nin();
    const int Nvol  = (*vk)[0].nvol();
    const int Nvol2 = Nvol / 2;

    for (int k = 0; k < 2 * m_Nsbt; ++k) {
      m_vk[k].reset(Nin, Nvol2, 1);
    }

    for (int k = 0; k < m_Nsbt; ++k) {
      m_ev[k] = (*ev)[k];
      m_index_eo.convertField(m_vk[k], (*vk)[k], 0);
      m_index_eo.convertField(m_vk[k + m_Nsbt], (*vk)[k], 1);
    }
  }

  //- setting parameter for low-mode subtraction
  m_prf.resize(m_Nsbt);
  for (int k = 0; k < m_Nsbt; ++k) {
    real_t sign_ev = m_ev[k] / fabs(m_ev[k]);
    m_prf[k] = (m_M0 - 0.5 * m_mq) * sign_ev
               - m_p0_parameter * m_ev[k];
  }

  const int Nsbt2 = 2 * m_Nsbt;

  m_u0c_e.resize(Nsbt2 * Nsbt2);
  m_u0cinv_e.resize(Nsbt2 * Nsbt2);
  m_u0c_o.resize(Nsbt2 * Nsbt2);
  m_u0cinv_o.resize(Nsbt2 * Nsbt2);

  if (m_Nsbt > 0) {
    Calc_Coeff_u0inv();
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::convert(AFIELD& v, const Field& w)
{
  if (!needs_convert()) {
    vout.crucial(m_vl, "%s: convert is not necessary.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier

  int Nex = w.nex();
  for (int ex = 0; ex < Nex; ++ex) {
    copy(m_v1, 0, w, ex);
    m_fopr_w->convert(m_v2, m_v1);
    copy(v, ex, m_v2, 0);
  }

#pragma omp barrier

}
//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::reverse(Field& v, const AFIELD& w)
{
  if (!needs_convert()) {
    vout.crucial(m_vl, "%s: reverse is not necessary.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp barrier

  int Nex = w.nex();
  for (int ex = 0; ex < Nex; ++ex) {
    copy(m_v1, 0, w, ex);
    m_fopr_w->reverse(m_v2, m_v1);
    copy(v, ex, m_v2, 0);
  }

#pragma omp barrier

}
//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::set_mode(const std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) m_mode = mode;

#pragma omp barrier

}
//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  int Nex = w.nex();
  assert(v.nex() == Nex);

  if(Nex == 1){
    m_fopr_w->mult_gm5(v, w);
  }else{
    for(int ex = 0; ex < Nex; ++ex){
      copy(m_v1, 0, w, ex);
      m_fopr_w->mult_gm5(m_v2, m_v1);
      copy(v, ex, m_v2, 0);
    }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D_eo") {
    DD_5d_eo(v, w, 1);
  } else if (m_mode == "DdagD_eo") {
    DdagD_eo(v, w);
  } else if (m_mode == "Ddag_eo") {
    DD_5d_eo(v, w, -1);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undefined\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D_eo") {
    DD_5d_eo(v, w, -1);
  } else if (m_mode == "DdagD_eo") {
    DdagD_eo(v, w);
  } else if (m_mode == "Ddag_eo") {
    DD_5d_eo(v, w, 1);
  } else {
    vout.crucial(m_vl, "Error at %s: mode undefined\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::DdagD_eo(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_Nin, m_Nvol2, m_Nex));
  assert(v.check_size(m_Nin, m_Nvol2, m_Nex));

  DD_5d_eo(m_t1, w, 1);
  DD_5d_eo(v, m_t1, -1);

}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::DD_5d_eo(AFIELD& v, const AFIELD& w,
                                        const int jd)
{
  assert(w.check_size(m_Nin, m_Nvol2, m_Nex));
  assert(v.check_size(m_Nin, m_Nvol2, m_Nex));

  if(jd == 1) {
    Mopr_5d_eo(m_t2, w, 1);
    LUprecond(v, m_t2, 1);
    Mopr_5d_eo(m_t2, v, 0);
    LUprecond(v, m_t2, 0);
    aypx(-1.0, v, w);    //    v = -v + w;
  }else if(jd == -1) {
    LUprecond(m_t2, w, 0);
    Mopr_5d_eo(v, m_t2, 1);
    LUprecond(m_t2, v, 1);
    Mopr_5d_eo(v, m_t2, 0);
    aypx(-1.0, v, w);    //    v = -v + w;
  }else{
    vout.crucial(m_vl, "Error at %s: illegal jd.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::Mopr_5d_eo(AFIELD& v, const AFIELD& w,
                                 const int ieo)
{
  //       ieo = 1: M_eo, 2: M_oe

  const int ieo1 = ieo;
  const int ieo2 = 1 - ieo;

  copy(m_w1, 0, w, 2 * m_Np);
  Proj_H_eo(ieo1, ieo2, m_z1, m_w1);

  m_z2.set(0.0);

  for (int j = 0; j < m_Np; ++j) {
    copy(m_w1, 0, w, 2 * j);
    if(ieo == 0){
      m_fopr_w->mult(m_v2, m_w1, "Deo");
    }else{
      m_fopr_w->mult(m_v2, m_w1, "Doe");
    }
    m_fopr_w->mult_gm5(m_v1, m_v2);
    copy(v, 2 * j, m_v1, 0);

    copy(m_w1, 0, w, 2 * j + 1);
    axpy(m_z2, m_p_sqrt[j], m_w1);

    if(ieo == 0){
      m_fopr_w->mult(m_v2, m_w1, "Deo");
    }else{
      m_fopr_w->mult(m_v2, m_w1, "Doe");
    }
    m_fopr_w->mult_gm5(m_v1, m_v2);
    scal(m_v1, -1.0);
    axpy(m_v1, m_p_sqrt[j], m_z1);
    copy(v, 2 * j + 1, m_v1, 0);
  }

  Proj_H_eo(ieo1, ieo2, m_v1, m_z2);
  copy(m_w1, 0, w, 2 * m_Np);
  Proj_L_mult_eo(ieo1, ieo2, m_z1, m_w1);
  axpy(m_v1, 1.0, m_z1);
  if(ieo == 0){
    m_fopr_w->mult(m_v2, m_w1, "Deo");
  }else{
    m_fopr_w->mult(m_v2, m_w1, "Doe");
  }
  m_fopr_w->mult_gm5(m_z1, m_v2);
  axpy(m_v1, m_p0_parameter, m_z1);
  copy(v, 2 * m_Np, m_v1, 0);
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::LUprecond(AFIELD& v, const AFIELD& w,
                                const int ieo)
{
  //      ieo=0: 1/M_ee, 1: 1/M_oo

  // --- L^-1 ---
  for (int ip = 0; ip < m_Np; ++ip) {
    int jx = 2 * ip;
    int jy = 2 * ip + 1;

    copy(m_z1, 0, w, jx);
    copy(m_z2, 0, w, jy);

    m_fopr_w->mult_gm5(m_v2, m_z1);
    axpy(m_z2, -m_rl[ip], m_v2);

    copy(m_t3, jx, m_z1, 0);
    copy(m_t3, jy, m_z2, 0);
  }

  const int j = 2 * m_Np;

  m_v1.set(0.0);
  for (int ip = 0; ip < m_Np; ++ip) {
    int jy = 2 * ip + 1;
    copy(m_z2, 0, m_t3, jy);
    axpy(m_v1, -m_sl[ip], m_z2);
  }
  m_fopr_w->mult_gm5(m_v2, m_v1);
  Proj_H_eo(ieo, ieo, m_v1, m_v2);
  copy(m_v2, 0, w, j);
  axpy(m_v1, 1.0, m_v2);
  copy(m_t3, j, m_v1, 0);

  // --- U^-1 ---
  mult_u0inv(m_v2, m_v1, ieo);
  copy(v, j, m_v2, 0);

  Proj_H_eo(ieo, ieo, m_v1, m_v2);
  copy(m_t3, j, m_v1, 0);

  for (int ip = 0; ip < m_Np; ++ip) {
    int jx = 2 * ip;
    int jy = 2 * ip + 1;

    copy(m_v2, 0, m_t3, jy);
    copy(m_v1, 0, m_t3, j);
    axpy(m_v2, -m_p_sqrt[ip], m_v1);
    m_fopr_w->mult_gm5(m_z2, m_v2);
    scal(m_z2, -1.0 / (1.0 + m_q_sqrt[ip] * m_q_sqrt[ip]));

    copy(v, jy, m_z2, 0);

    copy(m_v2, 0, m_t3, jx);
    axpy(m_v2, m_q_sqrt[ip], m_z2);
    m_fopr_w->mult_gm5(m_z1, m_v2);
    copy(v, jx, m_z1, 0);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::Proj_H_eo(const int ieo1, const int ieo2,
                                         AFIELD& v, const AFIELD& w)
{
  if (ieo1 == ieo2) {
    copy(v, w);
  } else {
    v.set(0.0);
  }

  for (int k = 0; k < m_Nsbt; ++k) {
    int j1 = k + ieo1 * m_Nsbt;
    int j2 = k + ieo2 * m_Nsbt;
    complex_t prod = dotc(m_vk[j2], w);
    axpy(v, -prod, m_vk[j1]);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::Proj_L_mult_eo(const int ieo1, const int ieo2,
                                              AFIELD& v, const AFIELD& w)
{
  v.set(0.0);

  for (int k = 0; k < m_Nsbt; ++k) {
    int j1 = k + ieo1 * m_Nsbt;
    int j2 = k + ieo2 * m_Nsbt;
    complex_t prod = dotc(m_vk[j2], w);
    prod *= m_prf[k];
    axpy(v, prod, m_vk[j1]);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::Calc_Coeff_u0inv()
{
  const int Nsbt  = m_Nsbt;
  const int Nsbt2 = 2 * m_Nsbt;

  const int Nin  = m_vk[0].nin();
  const int Nvol = m_vk[0].nvol();

  const real_t u0_a    = m_R_parameter + m_p0_parameter + m_u0;
  const real_t u0inv_a = 1.0 / u0_a;

  std::valarray<complex_t> c(Nsbt2 * Nsbt2);
  std::valarray<complex_t> cinv(Nsbt2 * Nsbt2);
  std::valarray<complex_t> vprod(Nsbt2 * Nsbt2);
  std::valarray<complex_t> W(Nsbt2 * Nsbt2);

  std::valarray<complex_t> c_src(Nsbt2);
  std::valarray<complex_t> c_x(Nsbt2);

  real_t a_r, a_i;

  for (int ieo = 0; ieo < 2; ++ieo) {
    //- Determine inner product of eigenvectors
    for (int i = 0; i < Nsbt; ++i) {
      for (int j = 0; j < i + 1; ++j) {
        int i2 = i + ieo * Nsbt;
        int j2 = j + ieo * Nsbt;

        complex_t a = dotc(m_vk[i2], m_vk[j2]);
        complex_t ac = cmplx(real(a), -imag(a));
	
        vprod[i + j * Nsbt2]                 = a;
        vprod[i + Nsbt + (j + Nsbt) * Nsbt2] = a;
        vprod[j + i * Nsbt2]                 = ac;
        vprod[j + Nsbt + (i + Nsbt) * Nsbt2] = ac;
      }
    }

    for (int i = 0; i < Nsbt; ++i) {
      for (int j = 0; j < i + 1; ++j) {
        int i2 = i + ieo * Nsbt;
        int j2 = j + ieo * Nsbt;

        m_fopr_w->mult_gm5(m_w1, m_vk[j2]);

        complex_t a = dotc(m_vk[i2], m_w1);
        complex_t ac = cmplx(real(a), -imag(a));

        vprod[i + (j + Nsbt) * Nsbt2] = a;
        vprod[i + Nsbt + j * Nsbt2]   = a;
        vprod[j + (i + Nsbt) * Nsbt2] = ac;
        vprod[j + Nsbt + i * Nsbt2]   = ac;
      }
    }

    //- definition of matrix c(i,j)
    c = cmplx(0.0, 0.0);
    for (int i = 0; i < Nsbt2; ++i) {
      c[i + i * Nsbt2] = cmplx(-m_u0, 0.0);
    }

    for (int i = 0; i < Nsbt; ++i) {
      real_t prf = (m_M0 - 0.5 * m_mq) * (m_ev[i] / fabs(m_ev[i]))
                   - m_p0_parameter * m_ev[i];
      c[i + (i + Nsbt) * Nsbt2] = cmplx(prf, 0.0);
    }

    for (int i = 0; i < Nsbt; ++i) {
      for (int j = Nsbt; j < Nsbt2; ++j) {
        c[i + j * Nsbt2] += cmplx(m_u0, 0.0) * vprod[i + j * Nsbt2];
      }
    }

    //- Definition of matrix W(i,j)
    W = cmplx(0.0, 0.0);
    for (int i = 0; i < Nsbt2; ++i) {
      W[i + i * Nsbt2] = cmplx(u0_a, 0.0);
    }

    for (int i = 0; i < Nsbt2; ++i) {
      for (int j = 0; j < Nsbt2; ++j) {
        for (int k = 0; k < Nsbt2; ++k) {
          W[i + j * Nsbt2] += c[i + k * Nsbt2] * vprod[k + j * Nsbt2];
        }
      }
    }

    //- Solve cinv
    for (int i = 0; i < Nsbt2; ++i) {
      for (int j = 0; j < Nsbt2; ++j) {
        c_src[j] = cmplx(-u0inv_a, 0.0) * c[j + i * Nsbt2];
      }

      Solv_Coeff_u0inv(c_x, W, c_src);

      for (int j = 0; j < Nsbt2; ++j) {
        cinv[j + i * Nsbt2] = c_x[j];
      }
    }

    if (ieo == 0) {
      for (int i = 0; i < Nsbt2 * Nsbt2; ++i) {
        m_u0c_e[i]    = c[i];
        m_u0cinv_e[i] = cinv[i];
      }
    } else {
      for (int i = 0; i < Nsbt2 * Nsbt2; ++i) {
        m_u0c_o[i]    = c[i];
        m_u0cinv_o[i] = cinv[i];
      }
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::mult_u0inv(AFIELD& v, const AFIELD& w,
                                 const int ieo)
{
  if (m_Nsbt == 0) {

    m_fopr_w->mult_gm5(v, w);
    scal(v, 1.0 / (m_R_parameter + m_p0_parameter + m_u0));

  } else {

    const int Nin   = w.nin();
    const int Nvol2 = w.nvol();

    const int Nsbt  = m_Nsbt;
    const int Nsbt2 = 2 * Nsbt;

    std::valarray<complex_t> prod_vb(Nsbt2), coeff(Nsbt2);
    std::valarray<complex_t> u0c(Nsbt2 * Nsbt2), u0cinv(Nsbt2 * Nsbt2);

    if (ieo == 0) {
      u0c    = m_u0c_e;
      u0cinv = m_u0cinv_e;
    } else {
      u0c    = m_u0c_o;
      u0cinv = m_u0cinv_o;
    }

    real_t a_r, a_i;
    for (int k = 0; k < Nsbt; ++k) {
      prod_vb[k] = dotc(m_vk[k + ieo * Nsbt], w);
    }
    m_fopr_w->mult_gm5(m_w1, w);
    for (int k = 0; k < Nsbt; ++k) {
      prod_vb[k + Nsbt] = dotc(m_vk[k + ieo * Nsbt], m_w1);
    }

    for (int i = 0; i < Nsbt2; ++i) {
      coeff[i] = cmplx(0.0, 0.0);
      for (int j = 0; j < Nsbt2; ++j) {
        coeff[i] += u0cinv[i + j * Nsbt2] * prod_vb[j];
      }
    }

    const real_t u0inv_a = 1.0 / (m_R_parameter + m_p0_parameter + m_u0);

    copy(m_w1, w);
    scal(m_w1, u0inv_a);

    for (int k = 0; k < Nsbt; ++k) {
      axpy(m_w1, coeff[k], m_vk[k + ieo * Nsbt]);
    }

    m_fopr_w->mult_gm5(v, m_w1);

    for (int k = 0; k < Nsbt; ++k) {
      axpy(v, coeff[k + Nsbt], m_vk[k + ieo * Nsbt]);
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::Solv_Coeff_u0inv(
                                    std::valarray<complex_t>& c_x,
                                    const std::valarray<complex_t>& W,
                                    const std::valarray<complex_t>& c_src)
{
  //- This is an implementation of CG solver.
  const int Nsbt  = m_Nsbt;
  const int Nsbt2 = 2 * m_Nsbt;

  assert(c_x.size() == Nsbt2);
  assert(c_src.size() == Nsbt2);
  assert(W.size() == (Nsbt2 * Nsbt2));

  const int    Niter = 100;
  const real_t Encg  = 1.e-32;

  const real_t ww    = norm_c(c_src);
  const real_t snorm = 1.0 / ww;

  std::valarray<complex_t> x(Nsbt2);
  std::valarray<complex_t> r(Nsbt2);
  std::valarray<complex_t> p(Nsbt2);
  std::valarray<complex_t> s(Nsbt2);
  std::valarray<complex_t> vt(Nsbt2);

  real_t rr, rr0;
  int    nconv = -1; // superficial initialization

  s = cmplx(0.0, 0.0);
  for (int i = 0; i < Nsbt2; ++i) {
    for (int j = 0; j < Nsbt2; ++j) {
      s[i] += conj(W[j + i * Nsbt2]) * c_src[j];
    }
  }
  x = s;
  r = s;
  mult_WdagW(s, W, x);

  r -= s;
  p  = r;
  rr = norm_c(r);

  vout.detailed(m_vl, "   init   %16.8e\n", rr * snorm);

  if (rr * snorm < Encg) goto converged;


  for (int iter = 0; iter < Niter; ++iter) {
    mult_WdagW(s, W, p);

    real_t pap   = innerprod_c(p, s);
    real_t alpha = rr / pap;

    x  += cmplx(alpha, 0.0) * p;
    r  -= cmplx(alpha, 0.0) * s;
    rr0 = rr;
    rr  = norm_c(r);

    vout.detailed(m_vl, "   %4d   %16.8e\n", iter, rr * snorm);

    if (rr * snorm < Encg) {
      nconv = iter;
      goto converged;
    }

    real_t beta = rr / rr0;

    p *= cmplx(beta, 0.0);
    p += r;
  }

  nconv = -1;

  vout.crucial(m_vl, "Error at %s: Not converged.\n", class_name.c_str());
  exit(EXIT_FAILURE);


converged:
  vout.detailed(m_vl, "  converged\n");

  s = cmplx(0.0, 0.0);
  for (int i = 0; i < Nsbt2; ++i) {
    for (int j = 0; j < Nsbt2; ++j) {
      s[i] += W[i + j * Nsbt2] * x[j];
    }
  }
  s -= c_src;

  real_t diff = norm_c(s);
  diff *= snorm;

  c_x = x;

  vout.general(m_vl, "  u0 solver: Nconv = %4d, diff = %12.4e\n", nconv, diff);
}


//====================================================================
template<typename AFIELD>
void AFopr_Overlap_5d<AFIELD>::mult_WdagW(
                                 std::valarray<complex_t>& v2,
                                 const std::valarray<complex_t>& W,
                                 const std::valarray<complex_t>& v1)
{
  const int size = v1.size();

  assert(v2.size() == size);
  assert(W.size() == (size * size));

  std::valarray<complex_t> vt(size);

  vt = cmplx(0.0, 0.0);
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      vt[i] += W[i + j * size] * v1[j];
    }
  }

  v2 = cmplx(0.0, 0.0);
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      v2[i] += conj(W[j + i * size]) * vt[j];
    }
  }
}


//====================================================================
template<typename AFIELD>
typename AFIELD::real_t AFopr_Overlap_5d<AFIELD>::norm_c(
                                     const std::valarray<complex_t>& v)
{
  const int size = v.size();

  real_t vv = 0.0;

  for (int i = 0; i < size; ++i) {
    vv += real(v[i]) * real(v[i]) + imag(v[i]) * imag(v[i]);
  }

  return vv;
}


//====================================================================
template<typename AFIELD>
typename AFIELD::real_t AFopr_Overlap_5d<AFIELD>::innerprod_c(
                                  const std::valarray<complex_t>& v,
                                  const std::valarray<complex_t>& w)
{
  const int size = v.size();

  real_t vw = 0.0;

  for (int i = 0; i < size; ++i) {
    vw += real(v[i]) * real(w[i]) + imag(v[i]) * imag(w[i]);
  }

  return vw;
}

//====================================================================
template<typename AFIELD>
double AFopr_Overlap_5d<AFIELD>::flop_count()
{
  //- Counting of floating point operations in giga unit.
  //  not implemented, yet.

  vout.general(m_vl, "Warning at %s: flop_count() has not been implemented.\n", class_name.c_str());

  const double gflop = 0.0;

  return gflop;
}


//============================================================END=====
