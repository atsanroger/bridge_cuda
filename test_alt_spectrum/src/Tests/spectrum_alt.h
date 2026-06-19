/*!
      @file    spectrum_alt.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 12:31:44 #$
      @version $LastChangedRevision: 2500 $
*/

#ifndef SPECTRUM_ALT_INCLUDED
#define SPECTRUM_ALT_INCLUDED

#include <vector>

#include  "lib/Parameters/commonParameters.h"
#include  "lib/Parameters/parameters.h"
#include  "lib/Field/field_G.h"
#include  "lib/IO/bridgeIO.h"
using Bridge::vout;

#include  "lib/Fopr/afopr.h"

//! Base class for spectrum measurement.

/*!
   This is a base class for spectrum measurements.
   Common operations to standard spectrum measurements,
   such as configuration setup, are provided.
                              [20 Apr 2017 H.Matsufuru]
 */

class Spectrum_alt
{
 public:
  static const std::string class_name;

 protected:
  Parameters params_all;
  Bridge::VerboseLevel m_vl;
  unique_ptr<Field_G> U;

 public:

  //! constructor
  Spectrum_alt() { init(); }

  //! destructor
  ~Spectrum_alt() { tidyup(); }

  //! configuration setup
  void setup_config(unique_ptr<Field_G>& U,
                    Parameters& params_all);

  //! gauge fixing
  void gauge_fixing(unique_ptr<Field_G>& U,
                    Parameters& params_all);

  //! gauge fixing as a test
  int gauge_fixing(std::string file_params, std::string run_mode);

  //  template<typename AFIELD>
  //  void mult_performance(AFopr<AFIELD>* fopr, int Nrepeat);

  //! check performance of linear algebra
  template<typename AFIELD>
  void check_linearAlgebra(int Nin, int Nvol, int Nex, int Nvec);

 private:

  //! initial setup
  void init();

  //! final tidy-up
  void tidyup();
};
#endif // SPECTRUM_ALT_INCLUDED
