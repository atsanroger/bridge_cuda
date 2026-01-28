/*!
        @file    force_F_Domainwall.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/

#ifndef FORCE_F_DOMAINWALL_INCLUDED
#define FORCE_F_DOMAINWALL_INCLUDED

#include "lib/Force/Fermion/force_F_Wilson_Nf2.h"

#include "lib/Fopr/fopr_Domainwall.h"
#include "lib/Field/index_lex.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for domain-wall fermions.

/*!
    At present, only the standard domain-wall setting is
    available.
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Domainwall : public Force
{
 public:
  static const std::string class_name;

 private:
  Index_lex *m_index;

  // parameters common to overlap fermion
  double m_mq;   // quark mass
  double m_M0;   // domain-wall height
  int    m_Ns;   // size of fifth-dimension
  std::vector<int> m_boundary;  //!< boundary conditions
  std::vector<double> m_b;  //!< coefficient parameter
  std::vector<double> m_c;  //!< coefficient parameter
  Bridge::VerboseLevel m_vl;   //!< verbose level

  std::string m_mode;

  Fopr_Domainwall *m_fopr_dw;
  Fopr_Wilson *m_fopr_w;
  Force_F_Wilson_Nf2 *m_force_w;

  //! initial setup.
  void init();

  //! final tidyup.
  void tidyup();

 public:

  //! constructor
  Force_F_Domainwall()
    : Force()
    { init(); }

  Force_F_Domainwall(const Parameters& params)
    : Force()
  {
    init();
    set_parameters(params);
  }

  //! destructor
  ~Force_F_Domainwall()
    { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double mq, const double M0, const int Ns,
                      const std::vector<int> bc,
                      const std::vector<double> b,
                      const std::vector<double> c);

  void set_config(Field *U);

  void set_mode(const std::string& mode) { m_mode = mode; }

  void force_udiv(Field& force, const Field& eta);

  void force_core1(Field& force, const Field& zeta, const Field& eta); // override default
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

  void force_udiv1_H(Field& force_, const Field& zeta,
                                    const Field& eta);

  void force_udiv1_Hdag(Field& force_, const Field& zeta,
                                       const Field& eta);

  void force_Himpl(Field& eta2, const Field& eta);

};
#endif
