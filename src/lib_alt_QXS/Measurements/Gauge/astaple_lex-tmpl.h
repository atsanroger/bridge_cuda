/*!
      @file    astaple_lex-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

template<typename AFIELD>
const std::string AStaple_lex<AFIELD>::class_name
                                          = "AStaple_lex<AFIELD>";
//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  int Nc = CommonParameters::Nc();
  m_Ndf = 2 * Nc * Nc;
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_shift = new ShiftAField_lex<AFIELD>(m_Ndf);

  m_Umu.reset(m_Ndf, m_Nvol, 1);

  m_v1.reset(m_Ndf, m_Nvol, 1);
  m_v2.reset(m_Ndf, m_Nvol, 1);
  m_v3.reset(m_Ndf, m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
	       class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::tidyup()
{
  delete m_shift;
}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::plaquette(real_t& plaq, const AFIELD& U)
{
  real_t plaqs, plaqt;
  plaq_s(plaqs, U);
  plaq_t(plaqt, U);

  plaq = 0.5 * (plaqs + plaqt);

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::plaquette(real_t& plaq, const Field& U)
{
  ThreadManager::assert_single_thread(class_name);

  int Nin  = U.nin();
  int Nvol = U.nvol();
  int Nex  = U.nex();
  AFIELD v1(Nin, Nvol, Nex);

  AIndex_lex<real_t,AFIELD::IMPL> index;
  convert_gauge(index, v1, U);

  real_t plaqs, plaqt;
  plaq_s(plaqs, v1);
  plaq_t(plaqt, v1);

  plaq = 0.5 * (plaqs + plaqt);

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::plaq_s(real_t& plaqs, const AFIELD& U)
{
  int nth = ThreadManager::get_num_threads();
  if(nth > 1){
    plaq_s_impl(plaqs, U);
  }else{
    plaq_s_omp(plaqs, U);
  }
}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::plaq_t(real_t& plaqt, const AFIELD& U)
{
  int nth = ThreadManager::get_num_threads();
  if(nth > 1){
    plaq_t_impl(plaqt, U);
  }else{
    plaq_t_omp(plaqt, U);
  }
}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::plaq_s_omp(real_t& plaqs, const AFIELD& U)
{

#pragma omp parallel
  {
    double plaq1;
    plaq_s_impl(plaq1, U);

#pragma omp master
    plaqs = plaq1;
  }

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::plaq_t_omp(real_t& plaqt, const AFIELD& U)
{

#pragma omp parallel
  {
    double plaq1;
    plaq_t_impl(plaq1, U);

#pragma omp master
    plaqt = plaq1;
  }

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::plaq_s_impl(real_t& plaqs, const AFIELD& U)
{
  int Nc       = CommonParameters::Nc();
  int Lvol     = CommonParameters::Lvol();

#pragma omp barrier

  real_t plaq = 0.0;

  for (int mu = 0; mu < m_Ndim-1; ++mu) {
    int nu = (mu + 1) % (m_Ndim-1);

    copy(m_Umu, 0, U, mu);
#pragma omp barrier

    upper(m_v3, U, mu, nu);
    //lower(m_v3, U, mu, nu);

    real_t plaq1 = dot(m_v3, m_Umu);
    plaq += plaq1;
  }

  plaqs = plaq/(Lvol * Nc * (m_Ndim-1));

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::plaq_t_impl(real_t& plaqt, const AFIELD& U)
{
  int Nc       = CommonParameters::Nc();
  int Lvol     = CommonParameters::Lvol();

#pragma omp barrier

  int mu = m_Ndim - 1;

  real_t plaq = 0.0;
  copy(m_Umu, 0, U, mu);
#pragma omp barrier

  for (int nu = 0; nu < m_Ndim-1; ++nu) {
    upper(m_v3, U, mu, nu);
    //lower(m_v3, U, mu, nu);

    real_t plaq1 = dot(m_v3, m_Umu);
    plaq += plaq1;
  }

  plaqt = plaq/(Lvol * Nc * (m_Ndim-1));

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::staple(AFIELD& w,
                                 const AFIELD& u, const int mu)
{
#pragma omp barrier

  w.set(0.0);

#pragma omp barrier

  for (int nu = 0; nu < m_Ndim; ++nu) {
    if (nu != mu) {
      upper(m_v3, u, mu, nu);
      axpy(w, real_t(1.0), m_v3);
#pragma omp barrier

      lower(m_v3, u, mu, nu);
      axpy(w, real_t(1.0), m_v3);
#pragma omp barrier
    }
  }

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::upper(AFIELD& c, const AFIELD& u,
                                const int mu, const int nu)
{
  // (1)  mu (2)
  //    +-->--+
  // nu |     |
  //   i+     +

  m_shift->backward(m_v1, 0, u, nu, mu);

  m_shift->backward(c, 0, u, mu, nu);

  QXS_Gauge::mult_Gnd(m_v2, 0, c, 0, m_v1, 0);

  QXS_Gauge::mult_Gnn(c, 0, u, nu, m_v2, 0);

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::lower(AFIELD& c, const AFIELD& u,
                                const int mu, const int nu)
{
  //    +     +
  // nu |     |
  //   i+-->--+
  //  (1)  mu (2)

  m_shift->backward(m_v2, 0, u, nu, mu);

  QXS_Gauge::mult_Gnn(m_v1, 0, u, mu, m_v2, 0);
  QXS_Gauge::mult_Gdn(m_v2, 0, u, nu, m_v1, 0);

  m_shift->forward(c, m_v2, nu);

}

//============================================================END=====
