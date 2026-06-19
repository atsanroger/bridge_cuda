/*!
        @file    staple_eo.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef STAPLE_EO_INCLUDED
#define STAPLE_EO_INCLUDED

#include "staple.h"

#include "Field/shiftField_eo.h"

#include "IO/bridgeIO.h"

//! Staple construction.

/*!
    Staple construction for even-odd gauge field.
                                    [28 Dec 2011 H.Matsufuru]
    Factory is introduced.          [24 Jan 2017 Y.Namekawa]
    Multi-threading is applied.
                                    [02 Dec 2021 H.Matsufuru]
 */

class Staple_eo : public Staple
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  std::string m_filename_output;

  int m_Nc, m_Nvol, m_Ndim;

  ShiftField_eo *m_shift;

  //  Field_G m_staple;
  Field_G m_v1, m_v2;
  Field_G m_w1, m_w2;
  Field_G m_Umu, m_Unu;

 public:

  Staple_eo() { init(); }

  Staple_eo(const Parameters& params) { init(params); }

  ~Staple_eo() { tidyup(); }

 public:
  void set_parameters(const Parameters& params);

  void get_parameters(Parameters& params) const;

  void upper(Field_G&, const Field_G&, const int, const int);

  void lower(Field_G&, const Field_G&, const int, const int);

  double plaq_s(const Field_G&);

  double plaq_t(const Field_G&);

  double plaquette(const Field_G&);

  void staple(Field_G&, const Field_G&, const int);

 private:
  // initial setup (obsolete).
  void init();

  // initial setup.
  void init(const Parameters& params);

  // final clean-up.
  void tidyup();

  double plaq_s_omp(const Field_G&);

  double plaq_t_omp(const Field_G&);

  double plaq_s_impl(const Field_G&);

  double plaq_t_impl(const Field_G&);


#ifdef USE_FACTORY
 private:
  static Staple *create_object()
  {
    return new Staple_eo();
  }

 public:
  static bool register_factory()
  {
    return Staple::Factory::Register("EvenOdd", create_object);
  }
#endif
};
#endif
