/*!
   @file    sample_HMC_quench.cpp
   @brief   Quenched (pure-gauge) HMC — generates a REAL, non-tiled thermalized
            SU(3) config. Stripped from sample_HMC.cpp: only Action_G, no fermion
            action / smearing / solver (so no MD-force blowups). CPU build.
*/
#include "bridge.h"
using namespace std;
using Bridge::vout;

namespace {
  const std::string test_name      = "sample_HMC_quench";
  const std::string parameter_file = "quench_HMC.yaml";
}

int quench_hmc(const Parameters& params_all)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const Parameters params_test     = params_all.lookup("Test_HMC_Wilson");
  const Parameters params_action_G = params_all.lookup("Action_G");
  const Parameters params_hmc      = params_all.lookup("HMC_Leapfrog");

  const string str_vlevel       = params_test.get_string("verbose_level");
  const Bridge::VerboseLevel vl  = vout.set_verbose_level(str_vlevel);

  const string str_gconf_status = params_test.get_string("gauge_config_status");
  const string str_gconf_read   = params_test.get_string("gauge_config_type_input");
  const string readfile         = params_test.get_string("config_filename_input");
  const string str_gconf_write  = params_test.get_string("gauge_config_type_output");
  const string writefile        = params_test.get_string("config_filename_output");
  const string str_rand_type    = params_test.get_string("random_number_type");
  const unsigned long seed      = params_test.get_unsigned_long("seed_for_random_number");
  int       i_conf              = params_test.get_int("trajectory_number");
  const int Ntraj               = params_test.get_int("trajectory_number_step");
  const int i_save_conf         = params_test.get_int("save_config_interval");

  const string str_action_G_type = params_action_G.get_string("action_type");

  RandomNumberManager::initialize(str_rand_type, seed);

  Field_G U(Nvol, Ndim);
  if (str_gconf_status == "Continue") {
    GaugeConfig(str_gconf_read).read(U, readfile);
  } else if (str_gconf_status == "Cold_start") {
    GaugeConfig("Unit").read(U);
  } else if (str_gconf_status == "Hot_start") {
    GaugeConfig("Random").read(U);
  } else {
    vout.crucial(vl, "Error at %s: bad gconf status %s\n", test_name.c_str(), str_gconf_status.c_str());
    exit(EXIT_FAILURE);
  }

  GaugeConfig gconf_write(str_gconf_write);

  unique_ptr<Action> action_G(Action::New(str_action_G_type, params_action_G));

  ActionList actions(1);
  actions.append(0, action_G.get());          // pure gauge — no fermion action
  std::vector<Director *> directors(0);        // no smearing director

  unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));
  HMC_Leapfrog hmc(actions, directors, rand.get(), params_hmc);

  Timer timer(test_name);
  timer.start();
  vout.general(vl, "Quenched HMC: Ntraj = %d\n", Ntraj);

  for (int traj = 0; traj < Ntraj; ++traj) {
    vout.general(vl, "\ntraj = %d\n", traj);
    hmc.update(U);
    if ((i_conf + traj + 1) % i_save_conf == 0) {
      std::string fn = FileUtils::generate_filename("%s-%06d", writefile.c_str(), (i_conf + traj + 1));
      gconf_write.write_file(U, fn);
    }
  }
  gconf_write.write_file(U, writefile);

  timer.report();
  RandomNumberManager::finalize();
  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  bridge_initialize(&argc, &argv);
  Parameters params = ParameterManager::read(parameter_file);
  bridge_setup(params.lookup("Main"));
  quench_hmc(params);
  bridge_finalize();
  return EXIT_SUCCESS;
}
