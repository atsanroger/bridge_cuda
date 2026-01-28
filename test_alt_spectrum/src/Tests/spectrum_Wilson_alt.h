/*!
        @file    spectrum_Wilson_alt.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-12-24 19:34:22 #$

        @version $LastChangedRevision: 2562 $
*/


#ifndef SPECTRUM_WILSON_ALT_INCLUDED
#define SPECTRUM_WILSON_ALT_INCLUDED

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
class Spectrum_Wilson_alt : public Spectrum_alt
{
 public:
  static const std::string class_name;

 private:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

 public:

  //! constructor
  Spectrum_Wilson_alt()
    : Spectrum_alt(), m_vl(CommonParameters::Vlevel())
  { init(); }

  //! destructor
  ~Spectrum_Wilson_alt() {}

  //
  //  int hadron_2ptFunction(std::string mode);

  //  template<Impl IMPL>
  int hadron_2ptFunction(std::string file_params, std::string mode);

  template<typename AFIELD>
  int check_operator(unique_ptr<Fopr>& fopr_ref,
                     Parameters& params_fopr);

  template<typename AFIELD>
  int check_operator_eo(unique_ptr<Fopr>& fopr_ref,
                        Parameters& params_fopr);

 private:

  //! initial setup
  void init();

  template<typename REALTYPE>
  void check(Fopr *, Source *, Parameters&);
};
#endif // SPECTRUM_WILSON_ALT_INCLUDED
