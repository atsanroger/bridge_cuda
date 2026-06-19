/*!
        @file    fermionFlow.cpp

        @brief

        @author  <Yusuke Taniguchi> tanigchi@het.ph.tsukuba.ac.jp

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#include "fermionFlow.h"

const std::string FermionFlow::class_name = "FermionFlow";

//====================================================================
void FermionFlow::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int              order_RK;
  double           Estep;
  int              Nprec;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_int("order_of_RungeKutta", order_RK);
  err += params.fetch_double("step_size", Estep);
  err += params.fetch_int("order_of_approx_for_exp_iP", Nprec);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "%s: fetch error, input parameter not found.\n", class_name.c_str());
    abort();
  }

  set_parameters(order_RK, Estep, Nprec, bc);
}


//====================================================================
void FermionFlow::get_parameters(Parameters& params) const
{
  params.set_int("order_of_RungeKutta", m_order_RK);
  params.set_double("step_size", m_Estep);
  params.set_int("order_of_approx_for_exp_iP", m_Nprec);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_int("number_of_steps", m_Nstep);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  return;
}


//====================================================================
void FermionFlow::set_parameters(const int order_RK,
                                 const double Estep, const int Nprec, const vector<int> bc)
{
  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  order_RK = %d\n", order_RK);
  vout.general(m_vl, "  Estep    = %10.6f\n", Estep);
  vout.general(m_vl, "  Nprec    = %d\n", Nprec);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(order_RK);
  err += ParameterCheck::square_non_zero(Estep);
  err += ParameterCheck::non_negative(Nprec);

  if (err) {
    vout.crucial(m_vl, "%s: parameter range check failed.\n", class_name.c_str());
    abort();
  }

  //- store values
  m_order_RK = order_RK; // order of Runge-Kutta
  m_Estep    = Estep;    // step size (SA)
  m_Nprec    = Nprec;    // order of approximation for e^{iP} (SA)

  int nsize = bc.size();
  assert(nsize == m_Ndim);
  m_boundary.resize(m_Ndim);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }
  m_gflow.set_parameters(m_order_RK, m_Estep, m_Nprec, 0, 1.0, 1.0);
}


//====================================================================
void FermionFlow::set_parameters(const int Nstep)
{
  //- print input parameters
  vout.general(m_vl, "  Nstep    = %d\n", Nstep);
  int err = 0;
  err += ParameterCheck::non_negative(Nstep);

  if (err) {
    vout.crucial(m_vl, "%s: parameter range check failed.\n", class_name.c_str());
    abort();
  }
  m_Nstep = Nstep;       // a number of steps (SA)
}


//====================================================================
void FermionFlow::laplacian(Field_F& v, const Field_F& f, Field_G& U)
{
  int     nvol = f.nvol();
  int     nex  = f.nex();
  Field_F trf(nvol, nex);
  Field_F trf2(nvol, nex);
  Field_F vt(nvol, 1);

  v.set(0.0);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_shift.backward(trf, f, m_boundary[mu], mu);
    for (int ex = 0; ex < f.nex(); ++ex) {
      mult_Field_Gn(vt, 0, U, mu, trf, ex);
      v.addpart_ex(ex, vt, 0);
      mult_Field_Gd(trf2, ex, U, mu, f, ex);
    }
    m_shift.forward(trf, trf2, m_boundary[mu], mu);
    axpy(v, 1.0, trf);
  }
  trf = f;
  scal(trf, 2.0 * m_Ndim);
  axpy(v, -1.0, trf);
}


//====================================================================
void FermionFlow::laplacian(Field_F_1spinor& v, const Field_F_1spinor& f, Field_G& U)
{
  int             nvol = f.nvol();
  int             nex  = f.nex();
  Field_F_1spinor trf(nvol, nex);
  Field_F_1spinor trf2(nvol, nex);
  Field_F_1spinor vt(nvol, 1);

  v.set(0.0);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_shift.backward(trf, f, m_boundary[mu], mu);
    for (int ex = 0; ex < nex; ++ex) {
      //      for (int isite = 0; isite < nvol; ++isite) {
      //        vt.set_vec(isite, 0, U.mat(isite, mu) * trf.vec(isite, ex));
      //      }
      mult_Field_Gn(vt, 0, U, mu, trf, ex);
      v.addpart_ex(ex, vt, 0);
      //      for (int isite = 0; isite < nvol; ++isite) {
      //        trf2.set_vec(isite, ex, U.mat_dag(isite, mu) * f.vec(isite, ex));
      //      }
      mult_Field_Gd(trf2, ex, U, mu, f, ex);
    }
    m_shift.forward(trf, trf2, m_boundary[mu], mu);
    axpy(v, 1.0, trf);
  }
  trf = f;
  scal(trf, 2.0 * m_Ndim);
  axpy(v, -1.0, trf);
}


//====================================================================
void FermionFlow::del_symmetric(int mu, Field_F& v, const Field_F& f, Field_G& U)
{
  int nvol = f.nvol();
  int nex  = f.nex();
  //Field_F w(nvol, nex);
  Field_F trf(nvol, nex);
  Field_F trf2(nvol, nex);

  assert(mu >= 0);
  assert(mu < m_Ndim);

  v.set(0.0);
  m_shift.backward(trf, f, m_boundary[mu], mu);
  for (int ex = 0; ex < nex; ++ex) {
    mult_Field_Gn(v, ex, U, mu, trf, ex);
    mult_Field_Gd(trf2, ex, U, mu, f, ex);
  }
  m_shift.forward(trf, trf2, m_boundary[mu], mu);
  axpy(v, -1.0, trf);
  scal(v, 0.5);
  //  v=(Field)w;
}


//====================================================================
void FermionFlow::del_symmetric(int mu, Field_F_1spinor& v, const Field_F_1spinor& f, Field_G& U)
{
  int             nvol = f.nvol();
  int             nex  = f.nex();
  Field_F_1spinor trf(nvol, nex);
  Field_F_1spinor trf2(nvol, nex);
  Field_F_1spinor vt(nvol, 1);

  v.set(0.0);
  m_shift.backward(trf, f, m_boundary[mu], mu);
  for (int ex = 0; ex < nex; ++ex) {
    for (int isite = 0; isite < nvol; ++isite) {
      vt.set_vec(isite, 0, U.mat(isite, mu) * trf.vec(isite, ex));
    }
    //      mult_Field_Gn(vt, 0, U, mu, trf, ex);
    v.addpart_ex(ex, vt, 0);
    for (int isite = 0; isite < nvol; ++isite) {
      trf2.set_vec(isite, ex, U.mat_dag(isite, mu) * f.vec(isite, ex));
    }
    //      mult_Field_Gd(trf, 0, U, mu, (Field_F)f, ex);
  }
  m_shift.forward(trf, trf2, m_boundary[mu], mu);
  axpy(v, -1.0, trf);
  scal(v, 0.5);
}


//====================================================================
void FermionFlow::evolve(Field_F& f, Field_G& U)
{
  assert(m_order_RK == 3);
  double  Eplaq = 0.0; // superficial initialization
  double  tt;
  Field_G U0((Field&)U);

  double t = 0.0;
  //- time evolution
  for (int i = 0; i < m_Nstep; ++i) {
    U = U0;
    m_action->set_config(&U); // give gauge conf pointer &U0 to d_action
    for (int j = 0; j < m_Nstep - i - 1; ++j) {
      m_gflow.evolve(t, U);
    }
    fermionFlow_3rd(f, U);
  }
}


//====================================================================
void FermionFlow::evolve(Field_F_1spinor& f, Field_G& U)
{
  assert(m_order_RK == 3);
  double  Eplaq = 0.0; // superficial initialization
  double  tt;
  Field_G U0((Field&)U);

  double t = 0.0;
  //- time evolution
  for (int i = 0; i < m_Nstep; ++i) {
    U = U0;
    m_action->set_config(&U); // give gauge conf pointer &U0 to d_action
    for (int j = 0; j < m_Nstep - i - 1; ++j) {
      m_gflow.evolve(t, U);
    }
    fermionFlow_3rd(f, U);
  }
}


//====================================================================
void FermionFlow::evolve(int mu, Field_F& f, Field_G& U)
{
  assert(m_order_RK == 3);
  double  t;
  Field_G U0((Field&)U);

  for (int i = 0; i < m_Nstep; ++i) {
    m_gflow.evolve(t, U);
  }
  Field_F fint(f);
  del_symmetric(mu, f, fint, U);

  //- time evolution
  for (int i = 0; i < m_Nstep; ++i) {
    U = U0;
    for (int j = 0; j < m_Nstep - i - 1; ++j) {
      m_gflow.evolve(t, U);
    }
    fermionFlow_3rd(f, U);
  }
}


//====================================================================
void FermionFlow::evolve2(Field_F_1spinor& f, const Field_G& U0, int measure, int gauge_store_interval)
{
  assert(m_order_RK == 3);
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_G U(Nvol, Ndim);

  int measurement_interval = m_Nstep / measure;

  double t = 0;
  //- time evolution
  for (int i = 0; i < measure; ++i) {
    for (int j = 0; j < measurement_interval; ++j) {
      int mm = (((measure - i) * measurement_interval - 1) / measurement_interval) / gauge_store_interval;
      t = mm * m_Estep;
      for (int mu = 0; mu < Ndim; ++mu) {
        U.setpart_ex(mu, U0, mu + Ndim * mm);
      }
      mm = ((measure - i) * measurement_interval - 1 - j) % (gauge_store_interval * measurement_interval);
      for (int k = 0; k < mm; ++k) {
        m_gflow.evolve(t, U);
      }
      fermionFlow_3rd(f, U);
    }
  }
}


//====================================================================
void FermionFlow::evolve_normal_order(Field_F_1spinor& f, Field_F& sf, Field_G& U)
{
  double Eplaq = 0.0;  // superficial initialization
  double tt;

  //- time evolution
  assert(m_order_RK == 3);
  for (int i = 0; i < m_Nstep; ++i) {
    fermionFlow_normal_order_3rd(f, sf, U);
  }
}


//====================================================================
void FermionFlow::fermionFlow_3rd(Field_F& f, Field_G& U)
{
  const double c3[] =
  { 0.25, -17.0 / 36.0, 8.0 / 9.0, 0.75, };
  Field_G W0((Field&)U);

  //- step 0
  m_action->force(m_iP, W0);                                   // calculate gradient of m_action Z_0
  mult_exp_Field_G(m_W_1, c3[0] * m_Estep, m_iP, W0, m_Nprec); // W_1=e^{Z_0/4}*U

  //- step 1
  scal(m_iP, c3[1]);             // -(17/36)Z_0
  m_action->force(m_iPP, m_W_1); // Z_1
  scal(m_iPP, c3[2]);            // (8/9)*Z_1
  axpy(m_iP, 1.0, m_iPP);        // (8/9)*Z_1-(17/36)*Z_0
  // W_2=e^{8*Z_1/9-17*Z_0/36}*W_1
  mult_exp_Field_G(m_U_0, m_Estep, m_iP, m_W_1, m_Nprec);

  Field_F lambda1(f.nvol(), f.nex());
  Field_F lambda2(f.nvol(), f.nex());
  //  Field lambda3=(Field&) f;
  Field_F scr(f.nvol(), f.nex());

  //  laplacian(lambda2, (Field&) f, m_U_0);
  laplacian(lambda2, f, m_U_0);
  scal(lambda2, m_Estep * c3[3]);

  laplacian(lambda1, lambda2, m_W_1);
  scal(lambda1, m_Estep * c3[2]);
  axpy(lambda1, 1.0, f);
  //  lambda1 += lambda3;

  scr = lambda2;
  scal(scr, -c3[2]);
  axpy(scr, 1.0, lambda1);
  laplacian(f, scr, W0);
  scal(f, m_Estep * c3[0]);
  axpy(f, 1.0, lambda1);
  axpy(f, 1.0, lambda2);
  //  vout.general(m_vl, "(%d %d %d %d) (%d) %.16f %.16f\n", 0, 0, 0, 0, 0, f.cmp_r(0, 0, 0, 0), f.cmp_i(0, 0, 0, 0));
  //  vout.general(m_vl, "(%d %d %d %d) (%d) %.16f %.16f\n", 0, 0, 0, 0, 0, W0.cmp_r(0, 0, 0), W0.cmp_i(0, 0, 0));
}


//====================================================================
void FermionFlow::fermionFlow_3rd(Field_F_1spinor& f, Field_G& U)
{
  const double c3[] =
  { 0.25, -17.0 / 36.0, 8.0 / 9.0, 0.75, };
  Field_G W0((Field&)U);

  //- step 0
  m_action->force(m_iP, W0);  // calculate gradient of m_action Z_0
  mult_exp_Field_G(m_W_1, c3[0] * m_Estep, m_iP, W0, m_Nprec);

  //- step 1
  scal(m_iP, c3[1]);             // -(17/36)Z_0
  m_action->force(m_iPP, m_W_1); // Z_1
  scal(m_iPP, c3[2]);            // (8/9)*Z_1
  axpy(m_iP, 1.0, m_iPP);        // (8/9)*Z_1-(17/36)*Z_0
  // W_2=e^{8*Z_1/9-17*Z_0/36}*W_1
  mult_exp_Field_G(m_U_0, m_Estep, m_iP, m_W_1, m_Nprec);

  Field_F_1spinor lambda1(f.nvol(), f.nex());
  Field_F_1spinor lambda2(f.nvol(), f.nex());
  Field_F_1spinor scr(f.nvol(), f.nex());

  laplacian(lambda2, f, m_U_0);
  scal(lambda2, m_Estep * c3[3]);

  laplacian(lambda1, lambda2, m_W_1);
  scal(lambda1, m_Estep * c3[2]);
  axpy(lambda1, 1.0, f);

  scr = lambda2;
  scal(scr, -c3[2]);
  axpy(scr, 1.0, lambda1);
  laplacian(f, scr, W0);
  scal(f, m_Estep * c3[0]);
  axpy(f, 1.0, lambda1);
  axpy(f, 1.0, lambda2);
}


//====================================================================
void FermionFlow::fermionFlow_normal_order_3rd(Field_F_1spinor& f, Field_F& sf, Field_G& U)
{
  const double c3[] =
  { 0.25, -17.0 / 36.0, 8.0 / 9.0, 0.75, -2.0 / 9.0, };

  Field_F_1spinor l0phi0(f.nvol(), f.nex());
  Field_F_1spinor phi1(f.nvol(), f.nex());
  Field_F_1spinor phi2(f.nvol(), f.nex());
  Field_F_1spinor scr(f.nvol(), f.nex());
  Field_F         l0sphi0(sf.nvol(), sf.nex());
  Field_F         sphi1(sf.nvol(), sf.nex());
  Field_F         sphi2(sf.nvol(), sf.nex());
  Field_F         sscr(sf.nvol(), sf.nex());

  //- step 1
  phi1  = f;
  sphi1 = (Field&)sf;
  laplacian(l0phi0, f, U);
  laplacian(l0sphi0, sf, U);
  axpy(phi1, m_Estep * c3[0], l0phi0);                    // phi1=phi0+1/4laplacian(0)phi0
  axpy(sphi1, m_Estep * c3[0], l0sphi0);                  // sphi1=sphi0+1/4laplacian(0)sphi0
  m_action->force(m_iP, U);                               // calculate gradient of m_action Z_0
  mult_exp_Field_G(U, c3[0] * m_Estep, m_iP, U, m_Nprec); // W_1=e^{Z_0/4}*U

  //- step 2
  phi2  = f;
  sphi2 = (Field&)sf;
  axpy(phi2, m_Estep * c3[4], l0phi0);
  axpy(sphi2, m_Estep * c3[4], l0sphi0);
  laplacian(scr, phi1, U);
  laplacian(sscr, sphi1, U);
  // phi2=phi0-2/9laplacian(0)phi0+8/9laplacian(1)phi1
  axpy(phi2, m_Estep * c3[2], scr);
  // sphi2=sphi0-2/9laplacian(0)sphi0+8/9laplacian(1)sphi1
  axpy(sphi2, m_Estep * c3[2], sscr);
  scal(m_iP, c3[1]);         // -(17/36)Z_0
  m_action->force(m_iPP, U); // Z_1
  scal(m_iPP, c3[2]);        // (8/9)*Z_1
  axpy(m_iP, 1.0, m_iPP);    // (8/9)*Z_1-(17/36)*Z_0
  // W_2=e^{8*Z_1/9-17*Z_0/36}*W_1
  mult_exp_Field_G(U, m_Estep, m_iP, U, m_Nprec);

  //- step 3
  laplacian(scr, phi2, U);
  laplacian(sscr, sphi2, U);
  axpy(phi1, m_Estep * c3[3], scr);   // phi3=phi1+3/4laplacian(2)phi2
  axpy(sphi1, m_Estep * c3[3], sscr); // sphi3=sphi1+3/4laplacian(2)sphi2
  f  = phi1;
  sf = sphi1;
  m_action->force(m_iPP, U); // Z_2
  scal(m_iPP, c3[3]);        // (3/4)*Z_2
  axpy(m_iPP, -1.0, m_iP);   // (3/4)*Z_2-(8/9)*Z_1+(17/36)*Z_0
  // V_out=e^{3*Z_2/4-8*Z_1/9+17*Z_0/36}*W_2
  mult_exp_Field_G(U, m_Estep, m_iPP, U, m_Nprec);

  //  vout.general(m_vl, "(%d %d %d %d) (%d) %.16f %.16f\n", 0, 0, 0, 0, 0, f.cmp_r(0, 0, 0, 0), f.cmp_i(0, 0, 0, 0));
  //  vout.general(m_vl, "(%d %d %d %d) (%d) %.16f %.16f\n", 0, 0, 0, 0, 0, W0.cmp_r(0, 0, 0), W0.cmp_i(0, 0, 0));
}


//============================================================END=====
