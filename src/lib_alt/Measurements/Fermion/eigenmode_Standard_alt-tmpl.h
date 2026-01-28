/*!
        @file    eigenmode_Standard_alt.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$
        @version $LastChangedRevision: 2422 $
*/

template<typename AFIELD>
const std::string Eigenmode_Standard_alt<AFIELD>::class_name
                                        = "Eigenmode_Standard_alt";
//====================================================================
template<typename AFIELD>
void Eigenmode_Standard_alt<AFIELD>::init(Parameters& params_fopr,
                                          Parameters& params_eigen)
{
  m_vl = CommonParameters::Vlevel();

  // this constructor assumes that the factories are available.
  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup.\n", class_name.c_str());

  string fopr_type = params_fopr.get_string("fermion_type");
  m_fopr = AFopr<AFIELD>::New(fopr_type, params_fopr);

  /*
  string eigensolver_type = params_eigensolver.get_string("eigensolver_type");
  m_eigensolver = AltEigensolver::New(eigensolver_type, m_fopr);
  m_eigensolver->set_parameters(params_eigensolver);
  */

  //  string eigensolver_type =
  //               params_eigensolver.get_string("eigensolver_type");
  m_eigensolver = new AEigensolver_IRLanczos<AFIELD, AFopr<AFIELD> >(m_fopr);
  m_eigensolver->set_parameters(params_eigen);


  m_Nk = params_eigen.get_int("number_of_wanted_eigenvectors");
  m_Np = params_eigen.get_int("number_of_working_eigenvectors");

  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void Eigenmode_Standard_alt<AFIELD>::tidyup()
{
  delete m_eigensolver;
  delete m_fopr;
}

//====================================================================
template<typename AFIELD>
void Eigenmode_Standard_alt<AFIELD>::set_config(Field *U)
{
  m_fopr->set_config(U);
}

//====================================================================
template<typename AFIELD>
void Eigenmode_Standard_alt<AFIELD>::calc_eigenmode()
// std::vector<double> TDa,
// std::vector<AField<double> >  vk, Nsbt, Nconv, b2)
{
  typedef typename AFIELD::real_t real_t;

  int Nm = m_Nk + m_Np;

  std::vector<real_t> TDa(Nm);
  std::vector<AFIELD>  vk(Nm);

  int NFin  = m_fopr->field_nin();
  int NFvol = m_fopr->field_nvol();
  int NFex  = m_fopr->field_nex();










}

//====================================================================
//============================================================END=====
