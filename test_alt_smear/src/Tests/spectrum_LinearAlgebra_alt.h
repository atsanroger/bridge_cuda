/*!
      @file    spectrum_LinearAlgebra_alt.h
      @brief
      @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-22 13:04:01 #$
      @version $LastChangedRevision: 2503 $
*/


#ifndef SPECTRUM_LINEARALGEBRA_ALT_INCLUDED
#define SPECTRUM_LINEARALGEBRA_ALT_INCLUDED

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
class Spectrum_LinearAlgebra_alt : public Spectrum_alt
{

public:
  static const std::string class_name;

private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

public:

  //! constructor
  Spectrum_LinearAlgebra_alt()
    : Spectrum_alt(), m_vl(CommonParameters::Vlevel())
    { init(); }

  //! destructor
  ~Spectrum_LinearAlgebra_alt(){}

  //
  int hadron_2ptFunction(std::string test_file, std::string mode);

 private:

  //! initial setup
  void init();

  template<typename REALTYPE>
    void check(Fopr*, Source*, Parameters&);

};
#endif // SPECTRUM_WILSON_ALT_INCLUDED
