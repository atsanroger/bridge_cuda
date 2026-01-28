/*!
        @file    hmc_Leapfrog.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef HMC_LEAPFROG_INCLUDED
#define HMC_LEAPFROG_INCLUDED

#include "action_list.h"
#include "integrator.h"
#include "langevin_Momentum.h"

#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC with single level leapfrog intetgrator.

/*!
    This class implements standartd HMC with a simple leapfrog
    molecular dynamics integrator.
    While more general integrator is now available, this class
    is easy to understand and convenient for a first test,
    and thus kept as it is.
                                     [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [03 Mar 2013 Y.Namekawa]
    Langevin step was moved to separate class.
                                     [07 May 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
*/


class HMC_Leapfrog
{
 public:
  static const std::string class_name;

 private:
  int m_Nmdc;
  int m_Nprec;
  bool m_Metropolis_test;
  double m_Estep;

  Langevin_Momentum *m_Langevin_P;
  Staple_lex *m_staple;

  std::vector<Action *> m_action;
  std::vector<Director *> m_director;
  RandomNumbers *m_rand;
  Bridge::VerboseLevel m_vl;

 public:

  //! constructor: with array of actions
  HMC_Leapfrog(std::vector<Action *> action,
               RandomNumbers *rand);

  //! constructor: with array of actions and directors
  HMC_Leapfrog(std::vector<Action *> action,
               std::vector<Director *> director,
               RandomNumbers *rand);

  //! constructor: with action_list
  HMC_Leapfrog(const ActionList& action_list,
               RandomNumbers *rand);

  //! constructor: with action_list and array of directors
  HMC_Leapfrog(const ActionList& action_list,
               std::vector<Director *> director,
               RandomNumbers *rand);

  //- constructor with paramsters
  //! constructor: with array of actions
  HMC_Leapfrog(std::vector<Action *> action,
               RandomNumbers *rand,
               const Parameters& params);

  //! constructor: with array of actions and directors
  HMC_Leapfrog(std::vector<Action *> action,
               std::vector<Director *> director,
               RandomNumbers *rand,
               const Parameters& params);

  //! constructor: with action_list
  HMC_Leapfrog(const ActionList& action_list,
               RandomNumbers *rand,
               const Parameters& params);

  //! constructor: with action_list and array of directors
  HMC_Leapfrog(const ActionList& action_list,
               std::vector<Director *> director,
               RandomNumbers *rand,
               const Parameters& params);

  //! destructor
  ~HMC_Leapfrog();

 private:
  // non-copyable
  HMC_Leapfrog(const HMC_Leapfrog&);
  HMC_Leapfrog& operator=(const HMC_Leapfrog&);

 public:
  void set_parameters(const Parameters& params);
  DEPRECATED
  void set_parameters(const double Estep, const int Nmdc, const int Nprec, const int Metropolis_test);  // backward compability
  void set_parameters(const double Estep, const int Nmdc, const int Nprec, const bool Metropolis_test);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  void get_parameters(Parameters& params) const;

  //  Metropolis_test: enabled if true.
  double update(Field_G&);

  double langevin(Field_G& iP, const Field_G& U);

  double calc_Hamiltonian(const Field_G& iP, const Field_G& U);
  double calcH_P(const Field_G& iP);

  void integrate(Field_G& iP, Field_G& U);

  void update_U(const double estep, const Field_G& iP, Field_G& U);
  void update_P(const double estep, Field_G& iP, const Field_G& U);
};
#endif
