/*!
      @file    asolver_SAP_MINRES-tmpl.h
      @brief   MinRes solver inside a SAP solver (Alt-version)
      @author  KANAMORI Issaku (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2#$
      @version $LastChangedRevision: 2499 $
 */

//====================================================================

template<typename AFIELD>
const std::string ASolver_SAP_MINRES<AFIELD>::class_name
  = "ASolver_SAP_MINRES";

//====================================================================
template<typename AFIELD>
void ASolver_SAP_MINRES<AFIELD>::init(void)
{
  ThreadManager::assert_single_thread(class_name);

  int nin  = m_fopr->field_nin();
  int nvol = m_fopr->field_nvol();
  int nex  = m_fopr->field_nex();

  m_r.reset(nin, nvol, nex);
  m_p.reset(nin, nvol, nex);

  m_p2_block.resize(m_block_index->coarse_nvol());
  m_alpha_block.resize(m_block_index->coarse_nvol());

#ifdef DEBUG_MINRES
  m_r2_block.resize(m_block_index->coarse_nvol());
#endif

  m_nconv = -1;
  calc_flop_each();
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_MINRES<AFIELD>::tidyup(void)
{
  // ThreadManager::assert_single_thread(class_name);
  // nothing is to be deleted.
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_MINRES<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Niter, Nrestart;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_int("maximum_number_of_restart", Nrestart);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int Niter2 = Niter * Nrestart;
  set_parameters(Niter2, Stop_cond);
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_MINRES<AFIELD>::set_parameters(const int Niter,
                                                const real_t Stop_cond)
{
  ThreadManager::assert_single_thread(class_name);

  m_Niter     = Niter;
  m_Stop_cond = Stop_cond;
  std::string prec = "double";
  if (sizeof(real_t) == 4) prec = "float";

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Precision: %s\n", prec.c_str());
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %16.8e\n", m_Stop_cond);
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_MINRES<AFIELD>::solve(AFIELD& x, const AFIELD& b,
                                       int& Nconv, real_t& diff, const int eo)
{
  assert(m_block_index);
  assert(m_fopr);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_block_index->coarse_nvol());

#pragma omp barrier // fine

  // init
  x.set(real_t(0.0));
  copy(m_r, b);

#pragma omp barrier // coarse

  for (int i = is; i < ns; ++i) {
    m_p2_block[i]    = real_t(1.0);
    m_alpha_block[i] = complex_t(0.0, 0.0);
  }

  // steps
  for (int iter = 0; iter < m_Niter; ++iter) {
    m_fopr->mult_sap(m_p, m_r, eo);

    block_norm2_eo(&(m_p2_block[0]), m_p, eo, *m_block_index);
    block_dotc_eo(&(m_alpha_block[0]), m_p, m_r, eo, *m_block_index);

    for (int i = is; i < ns; ++i) {
      if (m_p2_block[i] > 0.0) {
        m_alpha_block[i] /= m_p2_block[i];
      } else {
        // if the probe m_p is null, stop updating
        // as m_p is null, m_alpha is automatically zero.
        // m_alpha_block[i]=0.0;
      }
    }
    block_axpy_eo(x, &(m_alpha_block[0]), m_r, eo, real_t(1.0),
                  *m_block_index); // uses +m_alpha_block

    block_axpy_eo(m_r, &(m_alpha_block[0]), m_p, eo, real_t(-1.0),
                  *m_block_index); // used -m_alpha_block

#ifdef DEBUG_MINRES
    block_norm2_eo(&(m_r2_block[0]), m_r, eo, *m_block_index);

#pragma omp barrier

#pragma omp master
    {
      double r2 = 0.0;
      double p2 = 0.0;
      //      printf(" block  ie   r2   p2\n");
      for (int i = 0; i < m_block_index->coarse_nvol(); ++i) {
        //        printf(" %d    %d    %e   %e\n", i, m_block_index->block_eo(i),
        //                     m_r2_block[i], m_p2_block[i]);
        if (m_block_index->block_eo(i) == eo) {
          r2 += m_r2_block[i];
          p2 += m_p2_block[i];
        }
      }
      r2 = Communicator::reduce_sum(r2);
      p2 = Communicator::reduce_sum(p2);
      vout.general(m_vl, "iter=%d , sum(r2)=%e  sum(p2)=%e\n", iter, r2, p2);
    }// master
#pragma omp barrier
#endif // DEBUG_MINRES
  } // iter

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
double ASolver_SAP_MINRES<AFIELD>::flop_count()
{
  return m_Niter * m_flop_each;
}


//====================================================================
template<typename AFIELD>
void ASolver_SAP_MINRES<AFIELD>::calc_flop_each()
{
  int Nin  = m_fopr->field_nin();
  int Nvol = m_fopr->field_nvol();
  int Nex  = m_fopr->field_nex();
  int NPE  = CommonParameters::NPE();

  //  constexpr int Nc=3;
  //  constexpr int Nd=4;

  size_t nvol2     = m_block_index->coarse_nvol() / 2;
  size_t block_vol = m_block_index->block_nvol();

  /*

  int block_x=m_block_index->fine_lattice_size(0)/m_block_index->coarse_lattice_size(0);
  int block_y=m_block_index->fine_lattice_size(1)/m_block_index->coarse_lattice_size(1);
  int block_z=m_block_index->fine_lattice_size(2)/m_block_index->coarse_lattice_size(2);
  int block_t=m_block_index->fine_lattice_size(3)/m_block_index->coarse_lattice_size(3);

  constexpr int clover_site = 4*Nc*Nc*Nd*Nd;
  constexpr int hop_x_site= Nc*Nd*(4*Nc+2);
  constexpr int hop_y_site= Nc*Nd*(4*Nc+2);
  constexpr int hop_z_site= Nc*Nd*(4*Nc+2);
  constexpr int hop_t_site= Nc*Nd*(4*Nc+1);
  constexpr int accum_site= 4*Nc*Nd;

  // assumes Dirac rep.
  // todo move flop_fopr to the Clover (or Clover_SAP?)
  double flop_x=2.0*static_cast<double>(hop_x_site)*(block_x-1)*block_y*block_z*block_t*nvol2*NPE;
  double flop_y=2.0*static_cast<double>(hop_y_site)*block_x*(block_y-1)*block_z*block_t*nvol2*NPE;
  double flop_z=2.0*static_cast<double>(hop_z_site)*block_x*block_y*(block_z-1)*block_t*nvol2*NPE;
  double flop_t=2.0*static_cast<double>(hop_t_site)*block_x*block_y*block_z*(block_t-1)*nvol2*NPE;
  double flop_fopr
    =flop_x+flop_y+flop_z+flop_t+static_cast<double>(clover_site+accum_site)*nvol2*NPE;
  */

  double flop_fopr   = m_fopr->flop_count_sap();
  double flop_norm2  = 8.0 * Nin * block_vol * nvol2 * NPE;
  double flop_innerp = 8.0 * Nin * block_vol * nvol2 * NPE;
  double flop_axpy   = 8.0 * Nin * block_vol * nvol2 * NPE;

  double flop
    = (flop_fopr + flop_norm2 + flop_innerp + 2 * flop_axpy
       + 2 * nvol2 * NPE);

  m_flop_each = flop;
}


//============================================================END=====
