/*!
        @file    fprop.h

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FPROP_INCLUDED
#define FPROP_INCLUDED

#include "Fopr/fopr.h"

#include "IO/bridgeIO.h"

//! Base class for fermion propagator class family.

/*!
                                        [28 Dec 2011 H.Matsufuru]
    Introduce unique_ptr to avoid memory leaks.
                                        [21 Mar 2015 Y.Namekawa]
    Add flop_count.                     [ 8 Aug 2016 Y.Namekawa]
    m_mode is added as a protected member data.
    Methods added: invert(), set_mode(), reset_performance(),
      report_performance(), mult_performance().
                                       [22 Sep 2018 H.Matsufuru]
 */

class Fprop
{
 protected:
  Bridge::VerboseLevel m_vl;
  std::string m_mode;  // mode for inverter.

 public:
  Fprop()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Fprop() {}

 private:
  // non-copyable
  Fprop(const Fprop&);
  Fprop& operator=(const Fprop&);

 public:
  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void invert_D(Field&, const Field&, int&, double&)     = 0;
  virtual void invert_DdagD(Field&, const Field&, int&, double&) = 0;

  virtual void set_config(Field *) = 0;

  virtual double flop_count() = 0;

  //! invert accordingly to the mode. [22 Sep 2018 H.Matsufuru]
  virtual void invert(Field& x, const Field& b, int& nconv, double& diff)
  {
    if (m_mode == "D") {
      invert_D(x, b, nconv, diff);
    } else if (m_mode == "DdagD") {
      invert_DdagD(x, b, nconv, diff);
    } else {
      vout.crucial("Fprop: unsupported mode\n");
      exit(EXIT_FAILURE);
    }
  }

  //! set the mode for invert(). [22 Sep 2018 H.Matsufuru]
  virtual void set_mode(const std::string& mode) { m_mode = mode; }

  virtual void reset_performance()
  { }  // added [22 Sep 2018 H.Matsufuru]

  virtual void get_performance(double& flop_count, double& elapsed_time)
  { }  // added [07 Nov 2018 H.Matsufuru]

  virtual void report_performance()
  { }  // added [22 Sep 2018 H.Matsufuru]

  virtual void mult_performance(const std::string mode, const int Nrepeat)
  { }  // added [22 Sep 2018 H.Matsufuru]
};
#endif
