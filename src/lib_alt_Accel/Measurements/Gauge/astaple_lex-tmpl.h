/*!
      @file    astaple_lex-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-09-10 17:39:07 #$
      @version $LastChangedRevision: 2537 $
*/

template<typename AFIELD>
const std::string AStaple_lex<AFIELD>::class_name
                                          = "AStaple_lex<AFIELD>";
//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();

  int Nc = CommonParameters::Nc();
  m_Ndf = 2 * Nc * Nc;
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_shift = new ShiftAField_lex<AFIELD>(m_Ndf);

  m_Umu.reset(m_Ndf, m_Nvol, 1);

  m_v1.reset(m_Ndf, m_Nvol, 1);
  m_v2.reset(m_Ndf, m_Nvol, 1);
  m_v3.reset(m_Ndf, m_Nvol, 1);

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
void AStaple_lex<AFIELD>::plaq_s(real_t& plaqs, const AFIELD& U)
{
  int Nc       = CommonParameters::Nc();
  int Lvol     = CommonParameters::Lvol();

  real_t plaq = 0.0;

  for (int mu = 0; mu < m_Ndim-1; ++mu) {
    int nu = (mu + 1) % (m_Ndim-1);

    copy(m_Umu, 0, U, mu);

    upper(m_v3, U, mu, nu);

    real_t plaq1 = dot(m_v3, m_Umu);

    plaq += plaq1;
  }

  plaqs = plaq;
  plaqs = plaqs/(Lvol * Nc * (m_Ndim-1));

}


//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::plaq_t(real_t& plaqt, const AFIELD& U)
{
  int Nc       = CommonParameters::Nc();
  int Lvol     = CommonParameters::Lvol();

  int mu = m_Ndim - 1;

  real_t plaq = 0.0;
  copy(m_Umu, 0, U, mu);

  for (int nu = 0; nu < m_Ndim-1; ++nu) {
    upper(m_v3, U, mu, nu);
    real_t plaq1 = dot(m_v3, m_Umu);
    plaq += plaq1;
  }

  plaqt = plaq;
  plaqt = plaqt/(Lvol * Nc * (m_Ndim-1));

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::staple(AFIELD& W, const AFIELD& U,
                                                    const int mu)
{
#pragma omp barrier

  W.set(0.0);
#pragma omp barrier

  for (int nu = 0; nu < m_Ndim; ++nu) {
    if (nu != mu) {
      upper_th(m_v3, U, mu, nu);
      axpy(W, real_t(1.0), m_v3);
      lower_th(m_v3, U, mu, nu);
      axpy(W, real_t(1.0), m_v3);
    }
  }

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::upper(AFIELD& c, const AFIELD& U,
                                         const int mu, const int nu)
{
   upper_th(c, U, mu, nu);
}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::lower(AFIELD& c, const AFIELD& U,
                                         const int mu, const int nu)
{
   lower_th(c, U, mu, nu);
}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::upper_th(AFIELD& c, const AFIELD& U,
                                         const int mu, const int nu)
{
  // (1)  mu (2)
  //    +-->--+
  // nu |     |
  //   i+     +

  m_shift->backward(m_v1, 0, U, nu, mu);

  m_shift->backward(c, 0, U, mu, nu);

  mult_Gnd(m_v2, 0, c, 0, m_v1, 0);

  mult_Gnn(c, 0, U, nu, m_v2, 0);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AStaple_lex<AFIELD>::lower_th(AFIELD& c, const AFIELD& U,
                                         const int mu, const int nu)
{
  //    +     +
  // nu |     |
  //   i+-->--+
  //  (1)  mu (2)

  m_shift->backward(m_v2, 0, U, nu, mu);

  mult_Gnn(m_v1, 0, U, mu, m_v2, 0);
  mult_Gdn(m_v2, 0, U, nu, m_v1, 0);

  m_shift->forward(c, m_v2, nu);

#pragma omp barrier
}

//============================================================END=====
