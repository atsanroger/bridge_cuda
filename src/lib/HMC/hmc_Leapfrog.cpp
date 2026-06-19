/*!
        @file    hmc_Leapfrog.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "hmc_Leapfrog.h"

const std::string HMC_Leapfrog::class_name = "HMC_Leapfrog";

//====================================================================
//! constructor: with array of actions
HMC_Leapfrog::HMC_Leapfrog(std::vector<Action *> action,
                           RandomNumbers *rand)
  : m_vl(CommonParameters::Vlevel())
{
  m_action.resize(action.size());
  for (int i = 0; i < action.size(); ++i) {
    m_action[i] = action[i];
  }
  m_staple          = new Staple_lex;
  m_rand            = rand;
  m_Estep           = 0.0;
  m_Nmdc            = 0;
  m_Nprec           = 0;
  m_Metropolis_test = false;
  m_Langevin_P      = new Langevin_Momentum(m_rand);
}


//====================================================================
//! constructor: with array of actions and directors
HMC_Leapfrog::HMC_Leapfrog(std::vector<Action *> action,
                           std::vector<Director *> director,
                           RandomNumbers *rand)
  : m_vl(CommonParameters::Vlevel())
{
  m_action.resize(action.size());
  for (int i = 0; i < action.size(); ++i) {
    m_action[i] = action[i];
  }
  m_director.resize(director.size());
  for (int i = 0; i < director.size(); ++i) {
    m_director[i] = director[i];
  }
  m_staple          = new Staple_lex;
  m_rand            = rand;
  m_Estep           = 0.0;
  m_Nmdc            = 0;
  m_Nprec           = 0;
  m_Metropolis_test = false;
  m_Langevin_P      = new Langevin_Momentum(m_rand);
}


//====================================================================
//! constructor: with action_list
HMC_Leapfrog::HMC_Leapfrog(const ActionList& action_list,
                           RandomNumbers *rand)
  : m_vl(CommonParameters::Vlevel())
{
  ActionSet action_set = action_list.get_actions();

  m_action.resize(action_set.size());
  for (int i = 0; i < action_set.size(); ++i) {
    m_action[i] = action_set[i];
  }
  m_staple          = new Staple_lex;
  m_rand            = rand;
  m_Estep           = 0.0;
  m_Nmdc            = 0;
  m_Nprec           = 0;
  m_Metropolis_test = false;
  m_Langevin_P      = new Langevin_Momentum(m_rand);
}


//====================================================================
//! constructor: with action_list and array of directors
HMC_Leapfrog::HMC_Leapfrog(const ActionList& action_list,
                           std::vector<Director *> director,
                           RandomNumbers *rand)
  : m_vl(CommonParameters::Vlevel())
{
  ActionSet action_set = action_list.get_actions();

  m_action.resize(action_set.size());
  for (int i = 0; i < action_set.size(); ++i) {
    m_action[i] = action_set[i];
  }
  m_director.resize(director.size());
  for (int i = 0; i < director.size(); ++i) {
    m_director[i] = director[i];
  }
  m_staple          = new Staple_lex;
  m_rand            = rand;
  m_Estep           = 0.0;
  m_Nmdc            = 0;
  m_Nprec           = 0;
  m_Metropolis_test = false;
  m_Langevin_P      = new Langevin_Momentum(m_rand);
}


//====================================================================
//! constructor: with array of actions and parameters
HMC_Leapfrog::HMC_Leapfrog(std::vector<Action *> action,
                           RandomNumbers *rand,
                           const Parameters& params)
  : m_vl(CommonParameters::Vlevel())
{
  m_action.resize(action.size());
  for (int i = 0; i < action.size(); ++i) {
    m_action[i] = action[i];
  }
  m_staple          = new Staple_lex;
  m_rand            = rand;
  m_Estep           = 0.0;
  m_Nmdc            = 0;
  m_Nprec           = 0;
  m_Metropolis_test = false;
  m_Langevin_P      = new Langevin_Momentum(m_rand);

  set_parameters(params);
}


//====================================================================
//! constructor: with array of actions and directors and parameters
HMC_Leapfrog::HMC_Leapfrog(std::vector<Action *> action,
                           std::vector<Director *> director,
                           RandomNumbers *rand,
                           const Parameters& params)
  : m_vl(CommonParameters::Vlevel())
{
  m_action.resize(action.size());
  for (int i = 0; i < action.size(); ++i) {
    m_action[i] = action[i];
  }
  m_director.resize(director.size());
  for (int i = 0; i < director.size(); ++i) {
    m_director[i] = director[i];
  }
  m_staple          = new Staple_lex;
  m_rand            = rand;
  m_Estep           = 0.0;
  m_Nmdc            = 0;
  m_Nprec           = 0;
  m_Metropolis_test = false;
  m_Langevin_P      = new Langevin_Momentum(m_rand);

  set_parameters(params);
}


//====================================================================
//! constructor: with action_list and parameters
HMC_Leapfrog::HMC_Leapfrog(const ActionList& action_list,
                           RandomNumbers *rand,
                           const Parameters& params)
  : m_vl(CommonParameters::Vlevel())
{
  ActionSet action_set = action_list.get_actions();

  m_action.resize(action_set.size());
  for (int i = 0; i < action_set.size(); ++i) {
    m_action[i] = action_set[i];
  }
  m_staple          = new Staple_lex;
  m_rand            = rand;
  m_Estep           = 0.0;
  m_Nmdc            = 0;
  m_Nprec           = 0;
  m_Metropolis_test = false;
  m_Langevin_P      = new Langevin_Momentum(m_rand);

  set_parameters(params);
}


//====================================================================
//! constructor: with action_list and array of directors and parameters
HMC_Leapfrog::HMC_Leapfrog(const ActionList& action_list,
                           std::vector<Director *> director,
                           RandomNumbers *rand,
                           const Parameters& params)
  : m_vl(CommonParameters::Vlevel())
{
  ActionSet action_set = action_list.get_actions();

  m_action.resize(action_set.size());
  for (int i = 0; i < action_set.size(); ++i) {
    m_action[i] = action_set[i];
  }
  m_director.resize(director.size());
  for (int i = 0; i < director.size(); ++i) {
    m_director[i] = director[i];
  }
  m_staple          = new Staple_lex;
  m_rand            = rand;
  m_Estep           = 0.0;
  m_Nmdc            = 0;
  m_Nprec           = 0;
  m_Metropolis_test = false;
  m_Langevin_P      = new Langevin_Momentum(m_rand);

  set_parameters(params);
}


//====================================================================
//! destructor
HMC_Leapfrog::~HMC_Leapfrog()
{
  delete m_staple;
  delete m_Langevin_P;
}


//====================================================================
void HMC_Leapfrog::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double Estep;
  int    Nmdc, Nprec;
  bool   Metropolis_test;
  double traj_length;

  int err = 0;
  err += params.fetch_int("number_of_steps", Nmdc);
  err += params.fetch_int("order_of_exp_iP", Nprec);
  err += params.fetch_bool("Metropolis_test", Metropolis_test);

  //- step_size or trajectory_length: either, or both consistently
  int isset_ssize   = params.fetch_double("step_size", Estep);
  int isset_trajlen = params.fetch_double("trajectory_length", traj_length);

  if ((isset_ssize != 0) && (isset_trajlen != 0)) {
    // neither set. error.
    ++err;
  } else if ((isset_ssize == 0) && (isset_trajlen != 0)) {
    // using Estep.
  } else if ((isset_ssize != 0) && (isset_trajlen == 0)) {
    // using trajectory_length.
    Estep = traj_length / Nmdc;
  } else {
    // both set. check consistency.
    double diff = Estep * Nmdc - traj_length;

    if (fabs(diff) >= CommonParameters::epsilon_criterion()) {
      ++err;
    }
  }

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found or not appropriate.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Estep, Nmdc, Nprec, Metropolis_test);
}


//====================================================================
void HMC_Leapfrog::get_parameters(Parameters& params) const
{
  params.set_int("number_of_steps", m_Nmdc);
  params.set_int("order_of_exp_iP", m_Nprec);
  params.set_bool("Metropolis_test", m_Metropolis_test);

  params.set_double("step_size", m_Estep);
  params.set_double("trajectory_length", m_Estep * m_Nmdc);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void HMC_Leapfrog::set_parameters(const double Estep, const int Nmdc, const int Nprec, const int Metropolis_test)
{
  vout.crucial(m_vl, "%s: warning: integer variable for Metroplis_test is obsolete. use boolean parameter.\n", class_name.c_str());

  return set_parameters(Estep, Nmdc, Nprec, (Metropolis_test == 0 ? false : true));
}


//====================================================================
void HMC_Leapfrog::set_parameters(const double Estep, const int Nmdc, const int Nprec, const bool Metropolis_test)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Number of actions: %d\n", m_action.size());
  vout.general(m_vl, "  Estep           = %10.6f\n", Estep);
  vout.general(m_vl, "  Nmdc            = %d\n", Nmdc);
  vout.general(m_vl, "  Nprec           = %d\n", Nprec);
  vout.general(m_vl, "  Metropolis_test = %s\n", Metropolis_test ? "true" : "false");

  //- range check
  // NB. Estep,Nmdc,Nprec,Metropolis_test == 0 is allowed.

  //- store values
  m_Estep           = Estep;           // step size (SA)
  m_Nmdc            = Nmdc;            // a number of steps (SA)
  m_Nprec           = Nprec;           // order of approximation for e^{iP} (SA)
  m_Metropolis_test = Metropolis_test; // Metropolis test on/off (SA)
}


//====================================================================
double HMC_Leapfrog::update(Field_G& Uorg)
{
  const int Nc = CommonParameters::Nc();     //color (SA)

  Field_G U(Uorg);                           // set original gauge conf. to U (SA)

  for (int i = 0; i < m_action.size(); ++i) {
    m_action[i]->set_config(&U); // set gauge conf. for each action (SA)
  }

  const int Nin  = U.nin();  // 2x(complex) 3(color) x 3(color) part (SA)
  const int Nvol = U.nvol(); // local volume (SA)
  const int Nex  = U.nex();  // direction mu (SA)

  Field_G iP(Nvol, Nex);     // define momentum iP for gauge fields (SA)


  vout.general(m_vl, "HMC (single leapfrog) start.\n");

  //- Langevin step
  const double H_total0 = langevin(iP, U); // calculate initial hamiltonian (SA)
  // H=p^2/s + S_G + all fermion contributions (SA)
  vout.general(m_vl, "  H_total(init)  = %18.8f\n", H_total0);

  const double plaq0 = m_staple->plaquette(U); // calculate plaquette of old U (SA)
  vout.general(m_vl, "  plaq(init)     = %18.8f\n", plaq0);

  //- initial Hamiltonian
  //    H_total0 = calc_Hamiltonian(iP,U);

  //- molecular dynamical integration
  integrate(iP, U); // MD step (SA)

  //- trial Hamiltonian
  const double H_total1 = calc_Hamiltonian(iP, U); // calculate final hamiltonian (SA)
  vout.general(m_vl, "  H_total(trial) = %18.8f\n", H_total1);

  const double plaq1 = m_staple->plaquette(U); // calculate plaquette of new U (SA)
  vout.general(m_vl, "  plaq(trial)    = %18.8f\n", plaq1);


  //- Metropolis test
  vout.general(m_vl, "Metropolis test.\n");

  const double diff_H           = H_total1 - H_total0;  // Delta H (SA)
  const double exp_minus_diff_H = exp(-diff_H);         // e^{-Delta H} (SA)

  double rand = -1.0;
  if (m_Metropolis_test) {
    rand = m_rand->get(); //  random number 0 <= rand < 1 (SA)
  }
  vout.general(m_vl, "  H(diff)       = %10.6f\n", diff_H);
  vout.general(m_vl, "  exp(-diff_H)  = %10.6f\n", exp_minus_diff_H);
  vout.general(m_vl, "  Random number = %10.6f\n", rand);

  if (rand <= exp_minus_diff_H) { // accepted
    vout.general(m_vl, "  Accepted\n");
    Uorg = U;                     // set new U to gauge conf. (SA)
  } else {                        // rejected
    vout.general(m_vl, "  Rejected\n");
  }


  const double plaq_final = m_staple->plaquette(Uorg);
  vout.general(m_vl, "  plaq(final)    = %18.8f\n", plaq_final);

  return plaq_final;
}


//====================================================================
double HMC_Leapfrog::langevin(Field_G& iP, const Field_G& U)
{
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();
  const int NcA  = Nc * Nc - 1;                     // Nc^2-1 (SA)

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  vout.general(m_vl, "Langevin step.\n");

  // discard caches
  for (int i = 0; i < m_director.size(); ++i) {
    m_director[i]->notify_linkv(); // tell director that gauge filed is modified.(SA)
  }

  //- kinetic term
  const double H_iP = m_Langevin_P->set_iP(iP);
  //                         calculate hamiltonian of momenta iP. (SA)
  vout.general(m_vl, "  H_kin      = %18.8f\n", H_iP);
  vout.general(m_vl, "  H_kin/dof  = %18.8f\n", H_iP / NcA / Nvol / NPE / Ndim);


  double H_actions = 0.0;
  for (int i = 0; i < m_action.size(); ++i) {
    H_actions += m_action[i]->langevin(m_rand);
    // calculate contribution to hamiltonian from each action. (SA)
  }

  const double H_total = H_iP + H_actions; // total hamiltonian. (SA)
  return H_total;
}


//====================================================================
double HMC_Leapfrog::calc_Hamiltonian(const Field_G& iP, const Field_G& U)
{
  const int Nin  = U.nin(); // Local volume (SA)
  const int Nvol = U.nvol();
  const int Nex  = U.nex();

  const int Nc  = CommonParameters::Nc();
  const int Nd  = CommonParameters::Nd();
  const int NcA = Nc * Nc - 1;

  const int NPE = CommonParameters::NPE();

  vout.general(m_vl, "Hamiltonian calculation.\n");

  // kinetic term
  const double H_iP = calcH_P(iP); // calculate hamiltonian for iP (SA)
  vout.general(m_vl, "  H_kin      = %18.8f\n", H_iP);
  vout.general(m_vl, "  H_kin/dof  = %18.8f\n", H_iP / NcA / Nvol / NPE / Nex);

  double H_actions = 0.0;
  for (int i = 0; i < m_action.size(); ++i) {
    // calculate contribution to hamiltonian from each action (SA)
    H_actions += m_action[i]->calcH();
  }

  const double H_total = H_iP + H_actions;
  return H_total;
}


//====================================================================
double HMC_Leapfrog::calcH_P(const Field_G& iP)
{
  //calculate hamiltonian for iP: |iP|^2/2 (SA)
  const double hn   = iP.norm();
  const double H_iP = 0.5 * hn * hn;

  return H_iP;
}


//====================================================================
void HMC_Leapfrog::integrate(Field_G& iP, Field_G& U)
{
  vout.general(m_vl, "Integration.\n");

  const double estep  = m_Estep;
  const double estep2 = 0.5 * m_Estep;

  // Initial half step of update of h
  if (m_Nmdc > 0) {
    const int imdc = 0;
    vout.general(m_vl, "  imdc = %d\n", imdc);
    update_P(estep2, iP, U); // update momentum iP at first half step (SA)
  }

  // Molecular dynamics step
  for (int imdc = 1; imdc < m_Nmdc + 1; imdc++) {
    vout.general(m_vl, "  imdc = %d\n", imdc);

    update_U(estep, iP, U); // update gauge field U at full step (SA)

    if (imdc < m_Nmdc) {
      update_P(estep, iP, U);  // update momentum iP at full step (SA)
    } else {
      update_P(estep2, iP, U); // update momentum iP at last half step (SA)
    }
  }  // here imdc loop ends
}


//====================================================================
void HMC_Leapfrog::update_P(const double estep, Field_G& iP, const Field_G& U)
{
  const int Nin  = U.nin();
  const int Nvol = U.nvol();
  const int Nex  = U.nex();
  const int Nc   = CommonParameters::Nc();

  Field force(Nin, Nvol, Nex);

  force.set(0.0);

  Field force1(Nin, Nvol, Nex);
  // double H_actions = 0.0;
  for (int i = 0; i < m_action.size(); ++i) {
    m_action[i]->force(force1); // calculate force for each action (SA)
    axpy(force, estep, force1);
  }

  // iP = iP + step-size * force (SA)
  axpy(iP, 1.0, force);
}


//====================================================================
void HMC_Leapfrog::update_U(const double estep, const Field_G& iP, Field_G& U)
{
  const int Nvol = U.nvol();
  const int Nex  = U.nex();
  const int Nc   = CommonParameters::Nc();

  Mat_SU_N u0(Nc), u1(Nc), u2(Nc);
  Mat_SU_N h1(Nc);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      u0 = U.mat(site, ex);
      u1 = U.mat(site, ex);
      h1 = iP.mat(site, ex);

      for (int iprec = 0; iprec < m_Nprec; ++iprec) {
        double exf = estep / (m_Nprec - iprec); // step_size/(N-k) (SA)
        u2  = h1 * u1;                          // iP* u1 (SA)
        u2 *= exf;                              // step_size*iP*u1/(N-K)  (SA)
        u1  = u2;
        u1 += u0;                               // U + step_size*iP*u1/(N-K) (SA)
      }
      // u1 =sum_{k=0}^{N-1} (step_size * iP)^k/k!, N=m_Nprec (SA)
      u1.reunit();             // reunitarize u1 (SA)
      U.set_mat(site, ex, u1); // U = u1 (SA)
    }
  }

  for (int i = 0; i < m_director.size(); ++i) {
    m_director[i]->notify_linkv(); // notify all directors about update of U (SA)
  }
}


//====================================================================
//============================================================END=====
