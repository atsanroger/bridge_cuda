/*!
        @file    spectrum_Staggered_alt.h
        @brief
        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-22 13:04:01 #$
        @version $LastChangedRevision: 2503 $
*/


#ifndef SPECTRUM_STAGGERED_ALT_INCLUDED
#define SPECTRUM_STAGGERED_ALT_INCLUDED

#include <vector>

#include  "spectrum_alt.h"

#include  "lib_alt/alt_impl.h"

#include  "lib/Parameters/commonParameters.h"
#include  "lib/Parameters/parameters.h"
#include  "lib/Field/field_G.h"
#include  "lib/Fopr/fopr.h"
#include  "lib/IO/bridgeIO.h"
using Bridge::vout;

//class Fopr;
class Source;

template<Impl IMPL>
class Spectrum_Staggered_alt : public Spectrum_alt
{

public:
  static const std::string class_name;

private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  //unique_ptr<Field_G> U;

public:

  //! constructor
  Spectrum_Staggered_alt()
    : Spectrum_alt(), m_vl(CommonParameters::Vlevel())
    { init(); }

  //! destructor
  ~Spectrum_Staggered_alt(){}

  //
  int hadron_2ptFunction(std::string file_params,
                         std::string test_mode,
                         std::string run_mode)
  { hadron_2ptFunction_Cube(file_params, test_mode, run_mode); }

  int hadron_2ptFunction_Evenodd(std::string file_params,
                                 std::string test_mode,
                                 std::string run_mode);

  int hadron_2ptFunction_Cube(std::string file_params,
                              std::string test_mode,
                              std::string run_mode);

  //  int eigenspectrum(std::string mode);

  //  int eigenspectrum_alt(std::string mode);

 private:

  //! initial setup
  void init();

  int check_operator(Fopr*, Parameters&, Field_G *U);
  int check_operator_eo(Fopr*, Parameters&, Field_G *U);

};
#endif // SPECTRUM_WILSON_ALT_INCLUDED
