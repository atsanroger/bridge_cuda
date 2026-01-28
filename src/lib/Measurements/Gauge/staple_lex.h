/*!
        @file    staple_lex.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#ifndef STAPLE_LEX_INCLUDED
#define STAPLE_LEX_INCLUDED

#include "staple.h"

#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"

//! Staple construction.

/*!
    This class constructs the staple.
    While the originial version was written by J.Noaki,
    the present version is completely modified by H.Matsufuru
    except for the interface.
                                    [28 Dec 2011 H.Matsufuru]
    Thread-parallelized.
    void version of upper and lower functions added; these are
    faster than the versions returning Field_G object.
                                    [28 Sep 2013 H.Matsufuru]
    Add parameters for output.      [27 Jun 2016 Y.Namekawa]
    Factory is introduced.          [24 Jan 2017 Y.Namekawa]
 */

class Staple_lex : public Staple
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  std::string m_filename_output;

  ShiftField_lex *m_shift;

  Field_G m_Umu, m_Unu;     //!< working vectors
  Field_G m_w1, m_w2;       //!< working vectors
  Field_G m_v1, m_v2;       //!< working vectors

  int m_Nc, m_Nvol, m_Ndim;

 public:
  Staple_lex() { init(); }

  // optional
  Staple_lex(const Parameters& params) { init(params); }

  ~Staple_lex() { tidyup(); }

  //! setting parameters.
  void set_parameters(const Parameters& params);

  //! getting parameters.
  void get_parameters(Parameters& params) const;

  //! constructs upper staple in mu-nu plane.
  void upper(Field_G&, const Field_G&, const int mu, const int nu);

  //! constructs lower staple in mu-nu plane.
  void lower(Field_G&, const Field_G&, const int mu, const int nu);

  //! constructs staple in mu-direction (summing up nu-direction).
  void staple(Field_G&, const Field_G&, const int mu);

  //! calculates plaquette value.
  double plaquette(const Field_G&);

  //! calculates spatial plaquette value.
  double plaq_s(const Field_G&);

  //! calculates temporal plaquette value.
  double plaq_t(const Field_G&);


 private:

  void init(const Parameters& params);

  void init();

  void setup();

  void tidyup();

  double plaq_s_omp(const Field_G&);

  double plaq_t_omp(const Field_G&);

  double plaq_s_impl(const Field_G&);

  double plaq_t_impl(const Field_G&);


#ifdef USE_FACTORY
 private:
  static Staple *create_object()
  {
    return new Staple_lex();
  }

 public:
  static bool register_factory()
  {
    return Staple::Factory::Register("Lexical", create_object);
  }
#endif
};
#endif
