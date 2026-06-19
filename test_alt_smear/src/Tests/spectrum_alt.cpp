/*!
      @file    spectrum_alt.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-22 13:04:01 #$
      @version $LastChangedRevision: 2503 $
*/

#include "spectrum_alt.h"
#include "Tests/test.h"

#include "lib/Measurements/Gauge/gaugeFixing.h"
#include "lib/Measurements/Gauge/staple_lex.h"
#include "lib/IO/gaugeConfig.h"
#include "lib/Tools/randomNumberManager.h"

#include "Tests/job_Utils.h"

// class name
const std::string Spectrum_alt::class_name = "Spectrum_alt";
//====================================================================
void Spectrum_alt::init()
{
  m_vl = CommonParameters::Vlevel();

  //  vout.general(m_vl, "hello, spectrum.\n");

  //  RandomNumberManager::initialize("Mseries", 1234567UL);

}

//====================================================================
void Spectrum_alt::tidyup()
{

  //  RandomNumberManager::finalize();

}

//====================================================================
int Spectrum_alt::gauge_fixing(std::string file_params,
                               std::string run_mode)
{
  std::string test_name = "Gauge_fixing";

  vout.general("\n");
  vout.general("-------------------------------------------------"
               "-------------------\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("test file: %s\n", file_params.c_str());
  vout.general("run mode : %s\n", run_mode.c_str());

  string inputfile = "input";  // configuration numbers are stored

  int Nc   = CommonParameters::Nc();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(file_params);

  Parameters params_test = params_all.lookup("Test_Spectrum");

  string str_vlevel = params_test.get_string("verbose_level");
  Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

  string config_type_input  = params_test.get_string("gauge_config_type_input");
  string config_file_input = params_test.get_string("config_filename_input");

  vout.general(vl, "  config_type_input  = %s\n", config_type_input.c_str());
  vout.general(vl, "  config_file_input  = %s\n", config_file_input.c_str());

  string config_type_output = params_test.get_string("gauge_config_type_output");
  string config_file_output = params_test.get_string("config_filename_output");

  vout.general(vl, "  config_type_output = %s\n", config_type_output.c_str());
  vout.general(vl, "  config_file_output = %s\n", config_file_output.c_str());

  // setup random number manager
  RandomNumberManager::initialize("Mseries", 1234567UL);

  unique_ptr<GaugeConfig> gconf_read(new GaugeConfig(config_type_input));
  unique_ptr<GaugeConfig> gconf_write(new GaugeConfig(config_type_output));

  // Setup gauge configuration
  U.reset(new Field_G(Nvol, Ndim));

  Parameters params_gfix = params_all.lookup("GaugeFixing");

  unique_ptr<Staple> staple(new Staple_lex());

  unique_ptr<Timer> timer(new Timer(test_name));

  int iconf, nconf;
  if(run_mode == "test"){
    // iconf = 1;
    iconf = 200000;
    nconf = 1;
  }else{
    read_input(iconf, nconf, inputfile);
    vout.general(vl, "initial iconf   = %d\n", iconf);
    vout.general(vl, "number of confs = %d\n", nconf);
  }

  // ####  Execution main part  ####
  for(int i = 0; i < nconf; ++i){
    vout.general(vl, "\n");
    vout.general(vl, "iconf = %d\n", iconf);

    if(run_mode == "test"){
      gconf_read->read_file(*U, config_file_input);
    }else{
      gconf_read->read_file(*U, filename_config(config_file_input, iconf));
    }

    timer->start();

    gauge_fixing(U, params_gfix);

    if(run_mode == "job")
      gconf_write->write_file(*U, filename_config(config_file_output, iconf));

    ++iconf;
    if(run_mode == "job") write_input(iconf, nconf, inputfile);

    timer->report();

  }  // i-loop (configuration)


  RandomNumberManager::finalize();

  return 0;

}

//====================================================================
void Spectrum_alt::setup_config(unique_ptr<Field_G>& U,
                                Parameters& params_test)
{
  string gconf_status = params_test.get_string("gauge_config_status");
  string gconf_read   = params_test.get_string("gauge_config_type_input");
  string readfile     = params_test.get_string("config_filename_input");

  int err = 0;
  err += ParameterCheck::non_NULL(gconf_status);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameters have not been set\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Configuration setup:\n");
  vout.general(m_vl, "  gconf_status = %s\n", gconf_status.c_str());
  vout.general(m_vl, "  gconf_read   = %s\n", gconf_read.c_str());
  vout.general(m_vl, "  readfile     = %s\n", readfile.c_str());

#if defined USE_GROUP_SU3

  if (gconf_status == "Continue") {
    GaugeConfig(gconf_read).read(*U, readfile);
  } else if (gconf_status == "Cold_start") {
    GaugeConfig("Unit").read(*U);
  } else if (gconf_status == "Hot_start") {
    GaugeConfig("Random").read(*U);
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported gconf status \"%s\"\n",
                 class_name.c_str(), gconf_status.c_str());
    exit(EXIT_FAILURE);
  }

#else

  vout.general(m_vl, "  Nun-SU3 gauge group.\n");

  unique_ptr<FieldIO> gconf;
  if(gconf_read == "Binary"){
    gconf.reset(new FieldIO_Binary(IO_Format::Trivial));
  }else if(gconf_read == "Text"){
    gconf.reset(new FieldIO_Text(IO_Format::Trivial));
  }else if(gconf_read == "Text_4x4x4x8"){
    gconf.reset(new FieldIO_Text(IO_Format::Trivial));
  }else if(gconf_read == "NO_OUTPUT"){
    gconf.reset(new FieldIO_Text(IO_Format::Trivial));
  }else{
    vout.crucial(m_vl, "  read config type not supported.\n");
    exit(EXIT_FAILURE);
  }

  if (gconf_status == "Continue") {
    gconf->read_file(U.get(), readfile);
  } else if (gconf_status == "Cold_start") {
    U->set_unit();
  } else {
    vout.crucial(m_vl, "Error at Test_Spectrum: unsupported gconf"
                       " status \"%s\"\n", gconf_status.c_str());
    exit(EXIT_FAILURE);
  }

#endif

  Staple_lex staple;
  double plaq = staple.plaquette((Field_G)*U);
  vout.general(m_vl, "  plaq = %f\n", plaq);

  return;

}

//====================================================================
void Spectrum_alt::gauge_fixing(unique_ptr<Field_G>& U,
                                Parameters& params_gfix)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  //  Parameters params_gfix = params_all.lookup("GaugeFixing");
  string gfix_type = params_gfix.get_string("gauge_fixing_type");

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Gauge fixing:\n");
  vout.general(m_vl, "  Gauge fixing type = %s\n", gfix_type.c_str());

  if(gfix_type == "None") return;

  unique_ptr<Field_G> Ufix(new Field_G(Nvol, Ndim));
  unique_ptr<GaugeFixing> gfix(GaugeFixing::New(gfix_type));
  gfix->set_parameters(params_gfix);

  gfix->fix(*Ufix, *U);

  copy(*U, *Ufix);

  return;

}

//============================================================END=====
