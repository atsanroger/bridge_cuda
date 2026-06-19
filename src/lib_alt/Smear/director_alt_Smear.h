/*!
        @file    director_alt_Smear.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$
        @version $LastChangedRevision: 2529 $
*/

#ifndef DIRECTOR_ALT_SMEAR_INCLUDED
#define DIRECTOR_ALT_SMEAR_INCLUDED

#include <cassert>

#include "lib/Tools/director.h"
#include "lib/Smear/forceSmear.h"
#include "lib/Field/field_G.h"
#include "lib/Smear/asmear.h"
#include "lib/Smear/aprojection.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;


//! Manager of smeared configurations.

/*!
    Alt-code version of Director_Smear.
                                  [10 Feb 2023 H.Matsufuru]
 */

template<typename AFIELD>
class Director_alt_Smear : public Director
{

public:
  static const std::string class_name;
  typedef typename AFIELD::real_t real_t;

 private:

  Bridge::VerboseLevel m_vl;

  int m_Nsmear;                   //!< number of smearing to be applied
  Field_G *m_U;                   //!< original thin link var.
  std::vector<Field_G> m_Usmear;  //!< smeared configs.
  int m_status_linkv;             //!< set to zero when link var. is updated

  AProjection<AFIELD> *m_proj;    //!< projection operator.

  std::string m_smear_type;       //!< Smearing type
  ASmear<AFIELD>* m_smear;        //!< smearing operator

  AForceSmear<AFIELD>* m_force;

 public:

  //! constructor: only with Parameters object is available.
  Director_alt_Smear(const Parameters& params)
    : Director() { init(params); }

  //! reset parameters after construction.
  void set_parameters(const Parameters& params);

  //! reset number of smearing after construction.
  void set_parameters(const int Nsmear);

  //! get parameters
  void get_parameters(Parameters& params) const;

  //! get number of applied smearing operation
  int get_Nsmear() { return m_Nsmear; }

  //! get pointer to i-th smeared config (0th is original thin link)
  Field* getptr_smearedConfig(const int ismr);

  //! smeared config: to be discarded
  Field_G* get_config();

  //! intermediate config: to be discarded
  Field_G* get_config(const int ismr);

  //! intermediate config.
  void get_config(AFIELD& Usmr);

  //! intermediate config.
  void get_config(AFIELD& Usmr, const int ismr);

  //! set pointer to original thin link variable
  void set_config(Field *U);

  void force_udiv(Field_G&, const Field_G&, const Field_G&);

  //! to be called when configuration is updated
  void notify_linkv() {  m_status_linkv = 0; }

 private:

  //! initial setup with Parameters (new).
  void init(const Parameters& params);

  //! final clean-up.
  void tidyup();

  //! smearing is performed by calling a function of Smear object
  void smear();

};
#endif
