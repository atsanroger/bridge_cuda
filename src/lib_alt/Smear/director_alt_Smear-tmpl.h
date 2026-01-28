/*!
        @file    director_alt_Smear-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$
        @version $LastChangedRevision: 2529 $
*/

#include "lib_alt/Smear/director_alt_Smear.h"


template<typename AFIELD>
const std::string Director_alt_Smear<AFIELD>::class_name
                                     = "Director_alt_Smear<AFIELD>";

//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::init(const Parameters& params)

{
  m_vl = CommonParameters::Vlevel();
  string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  std::string m_smear_type  = params.get_string("smear_type");
  vout.general(m_vl, "smear_type   = %s\n", m_smear_type.c_str());

  if(m_smear_type != "None"){

    std::string proj_type = params.get_string("projection_type");
    vout.general(m_vl, "projection_type = %s\n", proj_type.c_str());
    m_proj = AProjection<AFIELD>::New(proj_type, params);

    m_smear = ASmear<AFIELD>::New(m_smear_type, m_proj, params);

  }else{
    m_proj = 0;
    m_smear = 0;
  }

  int err = 0;

  int Nsmear;
  err += params.fetch_int("number_of_smearing", Nsmear);

  if(err){
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Nsmear = Nsmear;

  vout.general(m_vl, "Nsmear = %d\n", m_Nsmear);
  std::string force_smear = params.get_string("force_smear");

  if(force_smear == "yes"){
    m_force = AForceSmear<AFIELD>::New(m_smear_type, m_proj, params);
    // m_force->set_parameters(params);
  }else{
    m_force = 0;
  }

  vout.general(m_vl, "verbose_level = %s\n", str_vlevel.c_str());

  const int Nvol = CommonParameters::Nvol();

  const int Ndim = CommonParameters::Ndim();

  m_Usmear.resize(m_Nsmear);
  if (m_Nsmear > 0) {
    for (int ismear = 0; ismear < m_Nsmear; ++ismear) {
      m_Usmear[ismear].reset(Nvol, Ndim);
    }
    vout.detailed(m_vl, " size of Usmear[ismear] was set.\n");
  }

  m_status_linkv = 0;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished\n", class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::tidyup()
{
  if(m_smear_type != "None"){
    delete m_smear;
    delete m_proj;
  }
}

//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int Nsmear;

  int err = 0;
  err += params.fetch_int("number_of_smearing", Nsmear);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Nsmear);
}


//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_string("smear_type", m_smear_type);

  if(m_smear_type != "None"){
    m_proj->get_parameters(params);
    m_smear->get_parameters(params);
  }

  params.set_int("number_of_smearing", m_Nsmear);
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::set_parameters(const int Nsmear)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  m_Nsmear = Nsmear;

  m_Usmear.resize(m_Nsmear);
  if (m_Nsmear > 0) {
    for (int ismr = 0; ismr < m_Nsmear; ++ismr) {
      m_Usmear[ismr].reset(Nvol, Ndim);
    }
    vout.paranoiac(m_vl, " size of Usmear[ismr] was set.\n");
  }

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  Nsmear = %d\n", m_Nsmear);

}


//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::set_config(Field *U)
{
  int ith = ThreadManager::get_thread_id();

  if(ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  if (m_status_linkv == 0) smear();
}


//====================================================================
template<typename AFIELD>
Field* Director_alt_Smear<AFIELD>::getptr_smearedConfig(const int ismr)
{
  assert(m_U != 0);

  if (ismr == 0) {
    return m_U;
  } else {
    return &m_Usmear[ismr - 1];
  }
}


//====================================================================
template<typename AFIELD>
Field_G* Director_alt_Smear<AFIELD>::get_config()
{
  if(m_Nsmear == 0){
    return m_U;
  }else{
    return &m_Usmear[m_Nsmear - 1];
  }

}

//====================================================================
template<typename AFIELD>
Field_G* Director_alt_Smear<AFIELD>::get_config(const int ismr)
{
  if(ismr == 0){
    return m_U;
  }else{
    return &m_Usmear[ismr-1];
  }

}

//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::get_config(AFIELD& Usmr)
{
#ifndef CORELIB_INSTANCE

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  if(m_Nsmear == 0){
    convert_gauge(index_lex, Usmr, *m_U);
  }else{
    convert_gauge(index_lex, Usmr, m_Usmear[m_Nsmear-1]);
  }

#endif

}

//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::get_config(AFIELD& Usmr, const int ismr)
{
#ifndef CORELIB_INSTANCE

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  if(m_Nsmear == 0){
    convert_gauge(index_lex, Usmr, *m_U);
  }else{
    convert_gauge(index_lex, Usmr, m_Usmear[ismr-1]);
  }

#endif
}

//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::force_udiv(Field_G& Sigma,
                                            const Field_G& Sigmap,
                                            const Field_G& U)
{
  m_force->force_udiv(Sigma, Sigmap, U);
}

//====================================================================
template<typename AFIELD>
void Director_alt_Smear<AFIELD>::smear()
{
  if (m_Nsmear > 0) {
    for (int ismr = 0; ismr < m_Nsmear; ++ismr) {
      if (ismr == 0) {
        m_smear->smear(m_Usmear[ismr], *m_U);
      } else {
        m_smear->smear(m_Usmear[ismr], m_Usmear[ismr - 1]);
      }
    }
  }

  ++m_status_linkv;
}

//============================================================END=====
