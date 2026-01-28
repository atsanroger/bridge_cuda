/*!
    @file    eigenvalue_alt.h
    @brief
    @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
             $LastChangedBy: matufuru $
    @date    $LastChangedDate:: 2022-12-16 15:57:38 #$
    @version $LastChangedRevision: 2422 $
*/

#ifndef EIGENVALUE_ALT_INCLUDED
#define EIGENVALUE_ALT_INCLUDED

#include <vector>
#include <string>

#include "spectrum_alt.h"

#include "lib_alt/alt_impl.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
//#include "lib/Field/field_F.h"
//#include "lib/Field/field_G.h"
//#include "lib/Fopr/fopr.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "job_Utils.h"

class Source;

template<typename AFIELD>
class Eigenvalue_alt : public Spectrum_alt
{

public:
  static const std::string class_name;

protected:
  using Spectrum_alt::m_vl;

public:

  //! constructor
  Eigenvalue_alt() : Spectrum_alt()
    { init(); }

  //! destructor
  ~Eigenvalue_alt(){}

  //! calculation of real eigenvalues.
  int measure(const std::string file_params,
              const std::string run_mode);

  //! calculation of real eigenvalues.
  int eigenvalue_Lanczos(const std::string file_params);

  //! calculation of complex eigenvalues.
  int eigenvalue_Arnoldi(const std::string file_params);

 private:

  //! initial setup
  void init();

};
#endif // EIGENVALUE_ALT_INCLUDED
