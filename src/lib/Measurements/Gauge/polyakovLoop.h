/*!
        @file    polyakovLoop.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef POLYAKOVLOOP_INCLUDED
#define POLYAKOVLOOP_INCLUDED

#include <cassert>
#include "Parameters/parameters.h"
#include "Field/field_G.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Polyakov loop measurement.

/*!
    This class determines the Polyakov loop for a given gauge
    configuration. Polyakov loop correlators are also planned to be
    measured, but still have not been implemented yet.
    For the latter case, set_parameters() is prepared (the definition
    of correlator constellation is the same as that of WilsonLoop).
                                       [22 Aug 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.               [14 Nov 2012 Y.Namekawa]
    Multithreaded.                     [29 Jul 2020 Y.Namekawa]
    Any number of threads is enabled.  [26 Dec 2021 H.Matsufuru]
 */


class PolyakovLoop
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  std::string m_filename_output;

  //! parameters set by user
  int m_Nspc_size;  //!< spatial size of loop
  int m_Ntype;      //!< number of measured loop-type

  //! working area
  Field_G m_P, m_Pcp1, m_Pcp2, m_Ut;

  //! internal data members (to be implemented).
  // int m_Ntype_max;  //!< maximum size of loop-type
  // int m_Nx_ext;     //!< size of extended gauge config.
  // int m_Ny_ext;     //!< size of extended gauge config.
  // int m_Nz_ext;     //!< size of extended gauge config.
  // int m_Nt_ext;     //!< size of extended gauge config.
  // int m_Nvol_ext;   //!< volume of extended gauge config.
  //
  // typedef std::vector<int> unitvec;
  // std::vector<unitvec> m_Nunit;
  // std::vector<int> m_Nmax;

 public:
  PolyakovLoop()
    : m_vl(CommonParameters::Vlevel()), m_Nspc_size(0), m_Ntype(0)
  {
    init();
  }

  PolyakovLoop(const Parameters& params)
    : m_vl(CommonParameters::Vlevel()), m_Nspc_size(0), m_Ntype(0)
  {
    init();
    set_parameters(params);
  }

  virtual ~PolyakovLoop() {}

 private:
  // non-copyable
  PolyakovLoop(const PolyakovLoop&);
  PolyakovLoop& operator=(const PolyakovLoop&);

 public:
  //! setting parameters: only for Polyakov loop correlators.
  virtual void set_parameters(const Parameters& params);
  void set_parameters(const int Nspc_size, const int Ntype);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void get_parameters(Parameters& params) const;

  //! Polyakov loop measurement
  dcomplex measure_ploop(const Field_G& U);

  //! Polyakov loop measurement (multi-threaded).
  void calc_ploop(dcomplex& ploop, const Field_G& U);

  //! Polyakov loop correlator measurement (to be implemented).
  // double measure_ploop_corr(const Field_G& U);

 private:

  //! initial setup independent of parameters.
  void init();
};
#endif
