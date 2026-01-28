/*!
      @file    spectrum_LinearAlgebra_alt-tmpl.h
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2022-12-16 15:57:38 #$
      @version $LastChangedRevision: 2422 $
*/

#include "spectrum_LinearAlgebra_alt.h"
#include "spectrum_alt-inc.h"

#include "Tests/test.h"

#include "lib/Field/field_F.h"
#include "lib/Field/field_G.h"

// class name
template<Impl IMPL>
const std::string Spectrum_LinearAlgebra_alt<IMPL>::class_name =
                                        "Spectrum_LinearAlgebra_alt";

//====================================================================
template<Impl IMPL>
void Spectrum_LinearAlgebra_alt<IMPL>::init()
{
  // do nothing.
}

//====================================================================
template<Impl IMPL>
int Spectrum_LinearAlgebra_alt<IMPL>::hadron_2ptFunction(
                             std::string test_file, std::string mode)
{
  const std::string test_name = class_name + "hadron_2ptFunction";

  vout.general("\n");
  vout.general("test name: %s\n", test_name.c_str());
  vout.general("test file: %s\n", test_file.c_str());

  // ####  parameter setup  ####
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  params_all = ParameterManager::read(test_file);

  Parameters params_test = params_all.lookup("Test_Spectrum");

  string str_vlevel = params_test.get_string("verbose_level");
  Bridge::VerboseLevel m_vl = vout.set_verbose_level(str_vlevel);
  vout.general(m_vl, "  vlevel       = %s\n", str_vlevel.c_str());

  // setup random number manager
  //  RandomNumberManager::initialize("Mseries", 1234567UL);

  // Setup gauge configuration 
  U.reset(new Field_G(Nvol, Ndim));
  setup_config(U, params_test);

  // Gauge fixing
  //  Parameters params_gfix = params_all.lookup("GaugeFixing");
  //  gauge_fixing(U, params_gfix);

  // fermion operator (reference)
  /*
  string str_fopr_type   = params_fopr.get_string("fermion_type");
  unique_ptr<Fopr> fopr(Fopr::New(str_fopr_type, str_gmset_type));
  fopr->set_parameters(params_fopr);
  fopr->set_config(U);
  */

  //qqqqq

  int Nin = 2 * Nc * Nd;
  int Nex = 1;
  int Nvec = 10;

  check_linearAlgebra<AField<double,IMPL> >(Nin, Nvol, Nex, Nvec);
  check_linearAlgebra<AField<float,IMPL> >(Nin, Nvol, Nex, Nvec);


  //  RandomNumberManager::finalize();

  /*
  if (do_check) {
    return Test::verify(result, expected_result);
  } else {
    vout.detailed(m_vl, "check skipped: expected_result not set.\n\n");
    return EXIT_SKIP;
  }
  */
  return EXIT_SKIP;

}

//============================================================END=====
