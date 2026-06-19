/*!
      @file    afopr_Domainwall_General_5din-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

template<typename AFIELD>
const std::string AFopr_Domainwall_General_5din<AFIELD>::class_name
                                     = "AFopr_Domainwall_General_5din";
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::init(const Parameters& params)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_repr = "Dirac";  // now only the Dirac repr is available.

  // switches
  //int req_comm = 1;  // set 1 if communication forced any time
  int req_comm = 0;  // set 1 if communication forced any time

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  m_Nvcd = 2 * Nc * Nd;
  m_Ndf  = 2 * Nc * Nc;

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nx = CommonParameters::Nx();
  m_Ny = CommonParameters::Ny();
  m_Nz = CommonParameters::Nz();
  m_Nt = CommonParameters::Nt();

  m_Nxv  = m_Nx/VLEN2;
  m_Nstv = m_Nvol/VLEN2;
  vout.general(m_vl, "  Nxv = %d  Nstv = %d\n", m_Nxv, m_Nstv);

  // setup verbose level
  string vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(vlevel);

  vout.general(m_vl, "Initialization of %s:\n", class_name.c_str());
  int err = 0;

  // setup kernel operator
  std::string kernel_type;
  err += params.fetch_string("kernel_type", kernel_type);
  if(err > 0){
    vout.crucial(m_vl, "Error at %s: kernel_type is not specified.\n",
		 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  Parameters params_kernel = params;
  double M0;
  err += params.fetch_double("domain_wall_height", M0);
  if(err > 0){
    vout.crucial(m_vl, "Error at %s: kernel_type is not specified.\n",
		 class_name.c_str());
    exit(EXIT_FAILURE);
  }
  m_M0 = real_t(M0);

  double kappa = 1.0/(8.0 - 2.0 * M0);
  params_kernel.set_double("hopping_parameter", kappa);

  // Factory is assumed to work.
  m_foprw = AFopr<AFIELD>::New(kernel_type, params_kernel);
  m_foprw->set_mode("D");

  m_Ns = 0;

  set_parameters(params);

  // gauge configuration.
  m_U.reset(m_Ndf, m_Nvol, m_Ndim);

  do_comm_any = 0;
  for(int mu = 0; mu < m_Ndim; ++mu){
    do_comm[mu] = 1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_Nbdsize.resize(m_Ndim);
  int Nvst = (m_Nvcd/2) * m_Ns * m_Nvol;
  m_Nbdsize[0] = Nvst / m_Nx;
  m_Nbdsize[1] = Nvst / m_Ny;
  m_Nbdsize[2] = Nvst / m_Nz;
  m_Nbdsize[3] = Nvst / m_Nt;

  setup_channels();

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::tidyup()
{
  delete m_foprw;
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for(int mu = 0; mu < m_Ndim; ++mu){

    size_t Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu] = Communicator::send_init(Nvsize, mu, -1);
    chsend_up[mu] = Communicator::send_init(Nvsize, mu,  1);
#ifdef USE_MPI
    chrecv_up[mu] = Communicator::recv_init(Nvsize, mu,  1);
    chrecv_dn[mu] = Communicator::recv_init(Nvsize, mu, -1);
#else
    void* buf_up = (void*)chsend_dn[mu]->ptr();
    chrecv_up[mu] = Communicator::recv_init(Nvsize, mu,  1, buf_up);
    void* buf_dn = (void*)chsend_up[mu]->ptr();
    chrecv_dn[mu] = Communicator::recv_init(Nvsize, mu, -1, buf_dn);
#endif

    if(do_comm[mu] == 1){
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::set_parameters(
                                             const Parameters& params)
{

  string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  string gmset_type;
  double mq, M0;
  int Ns;
  std::vector<int> bc;
  double b, c;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", gmset_type);

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double("coefficient_b", b);
  err += params.fetch_double("coefficient_c", c);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(mq, M0, Ns, bc, b, c);

  if(real_t(M0) != m_M0) set_kernel_parameters(params);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::set_parameters(
                                const double mq,
                                const double M0,
                                const int Ns,
                                const vector<int> bc,
                                const double b,
                                const double c)
{
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());

  m_M0 = real_t(M0);
  m_mq = real_t(mq);
  m_Ns = Ns;

  m_NinF = m_Nvcd * m_Ns;

  assert(bc.size() == m_Ndim);
  if(m_boundary.size() != m_Ndim) m_boundary.resize(m_Ndim);
  for(int mu = 0; mu < m_Ndim; ++mu){
    m_boundary[mu] = bc[mu];
  }

  if(m_b.size() != m_Ns){
    m_b.resize(m_Ns);
    m_c.resize(m_Ns);
  }
  for(int is = 0; is < m_Ns; ++is){
    m_b[is] = real_t(b);
    m_c[is] = real_t(c);
  }

  vout.general(m_vl, "  mq   = %8.4f\n",m_mq);
  vout.general(m_vl, "  M0   = %8.4f\n",m_M0);
  vout.general(m_vl, "  Ns   = %4d\n",m_Ns);
  for(int mu = 0; mu < m_Ndim; ++mu){
    vout.general(m_vl, "  boundary[%d] = %2d\n",mu,m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients b = %16.10f  c = %16.10f\n",
                                                  m_b[0], m_c[0]);

  set_precond_parameters();

  // working 5d vectors.
  if(m_w1.nin() != m_NinF){
    m_w1.reset(m_NinF, m_Nvol, 1);
    m_v1.reset(m_NinF, m_Nvol, 1);
    m_v2.reset(m_NinF, m_Nvol, 1);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::set_kernel_parameters(
                                             const Parameters& params)
{
  Parameters params_kernel = params;

  double M0;
  params.fetch_double("domain_wall_height", M0);

  double kappa = 1.0/(8.0 - 2.0 * M0);
  params_kernel.set_double("hopping_parameter", kappa);

  m_foprw->set_parameters(params_kernel);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::set_precond_parameters()
{

  if(m_dp.size() != m_Ns){
    m_dp.resize(m_Ns);
    m_dm.resize(m_Ns);
    m_e.resize(m_Ns-1);
    m_f.resize(m_Ns-1);
  }

  for(int is=0; is < m_Ns; ++is){
    m_dp[is] = 1.0 + m_b[is]*(4.0 - m_M0);
    m_dm[is] = 1.0 - m_c[is]*(4.0 - m_M0);
  }

  m_e[0] = m_mq * m_dm[m_Ns-1]/m_dp[0];
  m_f[0] = m_mq * m_dm[0];
  for(int is=1; is < m_Ns-1; ++is){
    m_e[is] = m_e[is-1]*m_dm[is-1]/m_dp[is];
    m_f[is] = m_f[is-1]*m_dm[is]/m_dp[is-1];
  }

  m_g = m_e[m_Ns-2]*m_dm[m_Ns-2];

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::set_coefficients(
                                  const std::vector<double> vec_b,
                                  const std::vector<double> vec_c)
{
  if(vec_b.size() != m_Ns || vec_c.size() != m_Ns){
    vout.crucial(m_vl, "%s: size of coefficient vectors incorrect.\n",
                 class_name.c_str());
  }

  vout.general(m_vl, "%s: coefficient vectors are set:\n",
                 class_name.c_str());

  for(int is = 0; is < m_Ns; ++is){
    m_b[is] = real_t(vec_b[is]);
    m_c[is] = real_t(vec_c[is]);
    vout.general(m_vl, "b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                       is, m_b[is], is, m_c[is]);
  }

  set_precond_parameters();

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::set_config(Field* u)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

#pragma omp parallel
 {
  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;
  convert_gauge(index_lex, m_U, *u);

#pragma omp barrier

  for(int mu = 0; mu < m_Ndim; ++mu){
    if(m_boundary[mu] != 1) set_boundary_config(m_U, mu);
  }

 }  // #pragma omp parallel

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::set_boundary_config(
                                                 AFIELD& U,
						 const int mu)
{
  int ipe[4], Nsize[4], Lsize[4];
  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;

  real_t bc = (real_t)(m_boundary[mu]);

  for(int i = 0; i < m_Ndim; ++i){
    ipe[i]   = Communicator::ipe(i);
    Nsize[i] = CommonParameters::Nsize(i);
    Lsize[i] = CommonParameters::Lsize(i);
  }

  int Nst = m_Nx * m_Ny * m_Nz * m_Nt;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nst);

  int j[4];
  for(int site = is; site < ns; ++site){
    int x = site % m_Nx;
    int yzt = site/m_Nx;
    int y = yzt % m_Ny;
    int z = (yzt/m_Ny) % m_Nz;
    int t = yzt/(m_Ny*m_Nz);
    j[0] = ipe[0] * Nsize[0] + x;
    j[1] = ipe[1] * Nsize[1] + y;
    j[2] = ipe[2] * Nsize[2] + z;
    j[3] = ipe[3] * Nsize[3] + t;

    if(j[mu] == Lsize[mu]-1){
      for(int in = 0; in < m_Ndf; ++in){
        int i = index_lex.idx_G(in, site, mu);
        real_t uv = bc * U.cmp(i);
        U.set(i, uv);
      }
    }

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::convert(AFIELD& v,
                                                  const Field& w)
{
  Index_lex_alt<real_t,AFIELD::IMPL> index;

  for(int site = 0; site < m_Nvol; ++site){
   for(int is = 0; is < m_Ns; ++is){
    for(int ivcd = 0; ivcd < m_Nvcd; ++ivcd){
      int in = ivcd + m_Nvcd * is;
      v.set(index.idx(in, m_NinF, site, 0),
                                   real_t(w.cmp(ivcd, site, is)));
    }
   }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::reverse(Field& v,
                                                  const AFIELD& w)
{
  Index_lex_alt<real_t,AFIELD::IMPL> index;

  for(int site = 0; site < m_Nvol; ++site){
   for(int is = 0; is < m_Ns; ++is){
    for(int ivcd = 0; ivcd < m_Nvcd; ++ivcd){
      int in = ivcd + m_Nvcd * is;
      v.set(ivcd, site, is,
                    double(w.cmp(index.idx(in, m_NinF, site, 0))) );
    }
   }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                                   std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  if(mode == "D"){
    D(v, w);
  }else if(mode == "Ddag"){
    Ddag(v, w);
  }else if(mode == "DdagD"){
    DdagD(v, w);
  }else if(mode == "D_prec"){
    D_prec(v, w);
  }else if(mode == "Ddag_prec"){
    Ddag_prec(v, w);
  }else if(mode == "DdagD_prec"){
    DdagD_prec(v, w);
  }else if(mode == "Prec"){
    Prec(v, w);
  }else if(mode == "Precdag"){
    Precdag(v, w);
  }else{
    vout.crucial(m_vl, "%s: undefined mode = %s\n",
                 class_name.c_str(), mode.c_str() );
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w,
                                                   std::string mode)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  if(mode=="D"){
    Ddag(v, w);
  }else if(mode=="Ddag"){
    D(v, w);
  }else if(mode=="DdagD"){
    DdagD(v, w);
  }else if(mode=="D_prec"){
    Ddag_prec(v, w);
  }else if(mode=="Ddag_prec"){
    D_prec(v, w);
  }else if(mode=="DdagD_prec"){
    DdagD_prec(v, w);
  }else if(mode=="Prec"){
    Precdag(v, w);
  }else if(mode=="Precdag"){
    Prec(v, w);
  }else{
    std::cout << "mode undeifined in AFopr_Domainwall_General_5din.\n";
    abort();
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::DdagD(AFIELD& v,
                                                  const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  D(m_v1, w);
  Ddag(v, m_v1);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::DdagD_prec(AFIELD& v,
                                                     const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nvol, m_Ns));
  assert(v.check_size(m_NinF, m_Nvol, m_Ns));

  D_prec(m_v1, w);
  Ddag_prec(v, m_v1);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::D_prec(AFIELD& v,
                                                   const AFIELD& w)
{
  #pragma omp barrier
  L_inv(v, w);
  U_inv(m_v2, v);
  D(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Ddag_prec(AFIELD& v,
                                                      const AFIELD& w)
{
  Ddag(v, w);
  Udag_inv(m_v2, v);
  Ldag_inv(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Prec(AFIELD& v,
                                                 const AFIELD& w)
{
  #pragma omp barrier
  L_inv(m_v2, w);
  U_inv(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Precdag(AFIELD& v,
                                                    const AFIELD& w)
{
  #pragma omp barrier
  Udag_inv(m_v2, w);
  Ldag_inv(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);  // working vector

  #pragma omp barrier

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    Vsimd_t vL[NCD], wL[NCD], yL[NCD], zL[NCD];

    for(int is = 0; is < m_Ns; ++is){

      int is_up = (is + 1) % m_Ns;
      real_t Fup = 0.5;
      if(is == m_Ns-1) Fup = -0.5 * m_mq;
      load_vec(wL, &wp2[Nin4*is_up], NCD);
      set_aPm5_dirac_vec(zL, Fup, wL, NC);

      int is_dn = (is - 1 + m_Ns) % m_Ns;
      real_t Fdn = 0.5;
      if(is == 0) Fdn = -0.5 * m_mq;
      load_vec(wL, &wp2[Nin4*is_dn], NCD);
      add_aPp5_dirac_vec(zL, Fdn, wL, NC);

      real_t FF1 = m_b[is] *(4.0 - m_M0) + 1.0;
      load_vec(wL, &wp2[Nin4*is], NCD);
      set_vec(vL, FF1, wL, NCD);
      real_t FF2 = m_c[is] *(4.0 - m_M0) - 1.0;
      axpy_vec(vL, FF2, zL, NCD);
      save_vec(&vp2[Nin4*is], vL, NCD);

      set_vec( yL, m_b[is], wL, NCD);
      axpy_vec(yL, m_c[is], zL, NCD);
      scal_vec(yL, real_t(-0.5), NCD);
      save_vec(&yp2[Nin4*is], yL, NCD);

    }

  }

  Dhop(vp, yp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);  // working vector

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  // v <- w  here gamma_5 is multiplied.
  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    Vsimd_t vL[NCD], yL[NCD];

    clear_vec(yL, NCD);
    for(int is = 0; is < m_Ns; ++is){
      load_mult_gm5_dirac_vec(vL, &wp2[Nin4*is], NC);
      save_vec(&vp2[Nin4*is], vL, NCD);
      save_vec(&yp2[Nin4*is], yL, NCD);
    }

  }

  Dhop(yp, vp);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];
    real_t *yp2 = &yp[Nin5 * site];

    Vsimd_t vL[NCD], yL[NCD], xL[NCD];

    real_t FF1, FF2;

    for(int is = 0; is < m_Ns; ++is){

      FF1 = m_b[is] *(4.0 - m_M0) + 1.0;

      load_mult_gm5_dirac_vec(xL, FF1, &wp2[Nin4*is], NC);
      load_vec(yL, &yp2[Nin4*is], NCD);
      axpy_vec(xL, real_t(-0.5)*m_b[is], yL, NCD);
      mult_gm5_dirac_vec(vL, xL, NC);

      int is_up = (is + 1) % m_Ns;
      FF2 = m_c[is_up] *(4.0 - m_M0) - 1.0;

      load_mult_gm5_dirac_vec(xL, FF2, &wp2[Nin4*is_up], NC);
      load_vec(yL, &yp2[Nin4*is_up], NCD);
      axpy_vec(xL, -real_t(0.5)*m_c[is_up], yL, NCD);

      real_t Fup = 0.5;
      if(is == m_Ns-1) Fup = -0.5 * m_mq;
      add_aPp5_dirac_vec(vL, Fup, xL, NC);

      int is_dn = (is - 1 + m_Ns) % m_Ns;
      FF2 = m_c[is_dn] *(4.0 - m_M0) - 1.0;

      load_mult_gm5_dirac_vec(xL, FF2, &wp2[Nin4*is_dn], NC);
      load_vec(yL, &yp2[Nin4*is_dn], NCD);
      axpy_vec(xL, -real_t(0.5)*m_c[is_dn], yL, NCD);

      real_t Fdn = 0.5;
      if(is == 0) Fdn = -0.5 * m_mq;
      add_aPm5_dirac_vec(vL, -Fdn, xL, NC);

      save_vec(&vp2[Nin4*is], vL, NCD);
    }

  }

  #pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Dhop(real_t* vp, real_t* wp)
{
  int Nin4  = VLEN * NC * ND;
  int Nin5  = Nin4 * m_Ns;
  int Nin4H = VLEN * NC * ND2;
  int Nin5H = Nin4H * m_Ns;
  int NvU   = m_Ndf * m_Nvol;

  int Nxy = m_Nxv * m_Ny;
  int Nxyz = m_Nxv * m_Ny * m_Nz;


  #pragma omp barrier

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  Dhop_1(vp, wp);

  #pragma omp barrier

  #pragma omp master
  {
    chset_send.start();
    chset_recv.start();
  }

  Dhop_b(vp, wp);

  #pragma omp master
  {
    chset_send.wait();
    chset_recv.wait();
  }

  #pragma omp barrier

  Dhop_2(vp, wp);

  #pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Dhop_1(real_t* vp, real_t* wp)
{
  int Nin4  = VLEN * NC * ND;
  int Nin5  = Nin4 * m_Ns;
  int Nin4H = VLEN * NC * ND2;
  int Nin5H = Nin4H * m_Ns;
  int NvU   = m_Ndf * m_Nvol;

  int Nxy = m_Nxv * m_Ny;
  int Nxyz = m_Nxv * m_Ny * m_Nz;


  #pragma omp barrier

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  real_t bufL[Nin5];


  for(int site = site0; site < site1; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;
    int ixy  = site % Nxy;
    int iy   = iyzt % m_Ny;
    int izt  = site/Nxy;
    int iz   = izt % m_Nz;
    int it   = izt / m_Nz;
    int ixyz = site % Nxyz;

    int idir;

    real_t *wp2 = &wp[Nin5 * site];

    idir = 0; // x-direction boundary
    if(do_comm[idir] == 1){

      if(ix == 0){
        real_t *buf1 = (real_t*)chsend_dn[idir]->ptr();
        int ibf = NVC * ND2 * m_Ns * iyzt;
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_xp1(&buf1[ibf + NVC*ND2*is], &wp2[Nin4*is], NC);
        }
      }

      if(ix == m_Nxv-1){
        real_t *buf1 = (real_t*)chsend_up[idir]->ptr();
        int ibf = NVC * ND2 * m_Ns * iyzt;
        real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_xm1(&buf1[ibf + NVC*ND2*is], up, &wp2[Nin4*is], NC);
        }
      }

    }

    idir = 1; // y-direction boundary
    if(do_comm[idir] == 1){

      if(iy == 0){
        real_t *buf1 = (real_t*)chsend_dn[idir]->ptr();
        int ibf = Nin5H * (ix + m_Nxv * izt);
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_yp1(&buf1[ibf+Nin4H*is], &wp2[Nin4*is], NC);
        }
      }

      if(iy == m_Ny-1){
        real_t *buf1 = (real_t*)chsend_up[idir]->ptr();
        int ibf = Nin5H * (ix + m_Nxv * izt);
        real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_ym1(&buf1[ibf+Nin4H*is], up, &wp2[Nin4*is], NC);
        }
      }

    }

    idir = 2; // z-direction boundary
    if(do_comm[idir] == 1){

      if(iz == 0){
        real_t *buf1 = (real_t*)chsend_dn[idir]->ptr();
        int ibf = Nin5H * (ixy + Nxy * it);
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_zp1(&buf1[ibf+Nin4H*is], &wp2[Nin4*is], NC);
        }
      }

      if(iz == m_Nz-1){
        real_t *buf1 = (real_t*)chsend_up[idir]->ptr();
        int ibf = Nin5H * (ixy + Nxy * it);
        real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_zm1(&buf1[ibf+Nin4H*is], up, &wp2[Nin4*is], NC);
        }
      }

    }

    idir = 3; // t-direction boundary
    if(do_comm[idir] == 1){

      if(it == 0){
        real_t *buf1 = (real_t*)chsend_dn[idir]->ptr();
        int ibf = Nin5H * ixyz;
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_tp1_dirac(&buf1[ibf+Nin4H*is], &wp2[Nin4*is], NC);
        }
      }

      if(it == m_Nt-1){
        real_t *buf1 = (real_t*)chsend_up[idir]->ptr();
        int ibf = Nin5H * ixyz;
        real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
        for(int is = 0; is < m_Ns; ++is){
          mult_wilson_tm1_dirac(&buf1[ibf+Nin4H*is], up, &wp2[Nin4*is], NC);
        }
      }

    }

  } // site loop (boundary)

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Dhop_b(real_t* vp, real_t* wp)
{
  int Nin4  = VLEN * NC * ND;
  int Nin5  = Nin4 * m_Ns;
  int Nin4H = VLEN * NC * ND2;
  int Nin5H = Nin4H * m_Ns;
  int NvU   = m_Ndf * m_Nvol;

  int Nxy = m_Nxv * m_Ny;
  int Nxyz = m_Nxv * m_Ny * m_Nz;


  #pragma omp barrier

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  real_t bufL[Nin5];

  for(int site = site0; site < site1; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;
    int ixy  = site % Nxy;
    int iy   = iyzt % m_Ny;
    int izt  = site/Nxy;
    int iz   = izt % m_Nz;
    int it   = izt / m_Nz;
    int ixyz = site % Nxyz;

    int idir, nei;

    real_t *wp2 = &wp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    real_t  z4[VLEN*NCD];
    Vsimd_t vL[NCD];

    idir = 0;

    if(ix < m_Nxv-1 || do_comm[idir] == 0){
      int ix2 = (ix + 1) % m_Nxv;
      nei = ix2 + m_Nxv * iyzt;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
      for(int is = 0; is < m_Ns; ++is){
	shift_vec2_bw(z4, &wp2[Nin4*is], &wpn[Nin4*is], NCD);
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_xpb(vL, up, z4, NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    if(ix > 0 || do_comm[idir] == 0){
      int ix2 = (ix - 1 + m_Nxv) % m_Nxv;
      nei = ix2 + m_Nxv * iyzt;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
      real_t *un = m_U.ptr(VLEN * NDF2 * nei  + NvU * idir);
      real_t uL[VLEN*NDF2];
      shift_vec2_fw(uL, up, un, NDF2);
      for(int is = 0; is < m_Ns; ++is){
        shift_vec2_fw(z4, &wp2[Nin4*is], &wpn[Nin4*is], NCD);
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_xmb(vL, uL, z4, NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    idir = 1; // y-direction bulk

    if(iy < m_Ny-1 || do_comm[idir] == 0){
      int iy2 = (iy + 1) % m_Ny;
      nei = ix + m_Nxv * (iy2 + m_Ny * izt);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_ypb(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    if(iy > 0 || do_comm[idir] == 0){
      int iy2 = (iy - 1 + m_Ny) % m_Ny;
      nei = ix + m_Nxv * (iy2 + m_Ny * izt);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_U.ptr(VLEN * NDF2 * nei + NvU * idir);

      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_ymb(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    idir = 2; // z-direction bulk

    if(iz < m_Nz-1 || do_comm[idir] == 0){
      int iz2 = (iz + 1) % m_Nz;
      nei = ixy + Nxy * (iz2 + m_Nz*it);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_zpb(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    if(iz > 0 || do_comm[idir] == 0){
      int iz2 = (iz - 1 + m_Nz) % m_Nz;
      nei = ixy + Nxy * (iz2 + m_Nz*it);
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_U.ptr(VLEN * NDF2 * nei + NvU * idir);
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_zmb(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    idir = 3; // t-direction bulk

    if(it < m_Nt-1 || do_comm[idir] == 0){
      int it2 = (it + 1) % m_Nt;
      nei = ixyz + Nxyz * it2;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_tpb_dirac(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

    if(it > 0 || do_comm[idir] == 0){
      int it2 = (it - 1 + m_Nt) % m_Nt;
      nei = ixyz + Nxyz * it2;
      real_t *wpn = &wp[Nin5 * nei];
      real_t *up = m_U.ptr(VLEN * NDF2 * nei + NvU * idir);
      for(int is = 0; is < m_Ns; ++is){
        load_vec(vL, &vp2[Nin4*is], NCD);
        mult_wilson_tmb_dirac(vL, up, &wpn[Nin4*is], NC);
        save_vec(&vp2[Nin4*is], vL, NCD);
      }
    }

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Dhop_2(real_t* vp, real_t* wp)
{
  int Nin4  = VLEN * NC * ND;
  int Nin5  = Nin4 * m_Ns;
  int Nin4H = VLEN * NC * ND2;
  int Nin5H = Nin4H * m_Ns;
  int NvU   = m_Ndf * m_Nvol;

  int Nxy = m_Nxv * m_Ny;
  int Nxyz = m_Nxv * m_Ny * m_Nz;


  #pragma omp barrier

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  real_t bufL[Nin5];

  for(int site = site0; site < site1; ++site){
    int ix   = site % m_Nxv;
    int iyzt = site/m_Nxv;
    int ixy  = site % Nxy;
    int iy   = iyzt % m_Ny;
    int izt  = site/Nxy;
    int iz   = izt % m_Nz;
    int it   = izt / m_Nz;
    int ixyz = site % Nxyz;

    int idir, nei;

    real_t *wp2 = &wp[Nin5 * site];
    real_t *vp2 = &vp[Nin5 * site];

    real_t  z4[VLEN*NCD];
    Vsimd_t vL[NCD];

    idir = 0;
    if(do_comm[idir] == 1){

      if(ix == m_Nxv-1){
        real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
        real_t *buf2 = (real_t*)chrecv_up[idir]->ptr();
        int ibf = NVC * ND2 * m_Ns * iyzt;
        for(int is = 0; is < m_Ns; ++is){
          shift_vec0_bw(z4, &wp2[Nin4*is], NCD);
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_xpb(vL, up, z4, NC);
          mult_wilson_xp2(vL, up, &buf2[ibf + NVC*ND2*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

      if(ix == 0){
        real_t *up = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
        real_t *buf2 = (real_t*)chrecv_dn[idir]->ptr();
        int ibf = NVC * ND2 * m_Ns * iyzt;
        real_t uL[VLEN*NDF2];
        shift_vec0_fw(uL, up, NDF2);
        for(int is = 0; is < m_Ns; ++is){
          shift_vec0_fw(z4, &wp2[Nin4*is], NCD);
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_xmb(vL, uL, z4, NC);
          mult_wilson_xm2(vL, &buf2[ibf + NVC*ND2*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

    }

    idir = 1; // y-direction boundary
    if(do_comm[idir] == 1){

      if(iy == m_Ny-1){
        real_t *buf2 = (real_t*)chrecv_up[idir]->ptr();
        int ibf = Nin5H * (ix + m_Nxv * izt);
        real_t *uL = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_yp2(vL, uL, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

      if(iy == 0){
        real_t *buf2 = (real_t*)chrecv_dn[idir]->ptr();
        int ibf = Nin5H * (ix + m_Nxv * izt);
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_ym2(vL, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

    }

    idir = 2; // y-direction boundary
    if(do_comm[idir] == 1){

      if(iz == m_Nz-1){
        real_t *buf2 = (real_t*)chrecv_up[idir]->ptr();
        int ibf = Nin5H * (ixy + Nxy * it);
        real_t *uL = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_zp2(vL, uL, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

      if(iz == 0){
        real_t *buf2 = (real_t*)chrecv_dn[idir]->ptr();
        int ibf = Nin5H * (ixy + Nxy * it);
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_zm2(vL, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

    }

    idir = 3; // y-direction boundary
    if(do_comm[idir] == 1){

      if(it == m_Nt-1){
        real_t *buf2 = (real_t*)chrecv_up[idir]->ptr();
        int ibf = Nin5H * ixyz;
        real_t *uL = m_U.ptr(VLEN * NDF2 * site + NvU * idir);
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_tp2_dirac(vL, uL, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

      if(it == 0){
        real_t *buf2 = (real_t*)chrecv_dn[idir]->ptr();
        int ibf = Nin5H * ixyz;
        for(int is = 0; is < m_Ns; ++is){
          load_vec(vL, &vp2[Nin4*is], NCD);
          mult_wilson_tm2_dirac(vL, &buf2[ibf+Nin4H*is], NC);
          save_vec(&vp2[Nin4*is], vL, NCD);
        }
      }

    }

  } // site loop (boundary)


  #pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::L_inv(AFIELD& v,
                                                  const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t vL[NCD], xL[NCD], yL[NCD];

    load_vec(vL, &wp2[0], NCD);
    save_vec(&vp2[0], vL, NCD);
    set_vec(yL, m_e[0], vL, NCD);

    for(int is = 1; is < m_Ns-1; ++is){
      copy_vec(xL, vL, NCD);
      load_vec(vL, &wp2[Nin4*is], NCD);
      add_aPp5_dirac_vec(vL, real_t(0.5)*m_dm[is]/m_dp[is-1], xL, NC);
      save_vec(&vp2[Nin4*is], vL, NCD);
      axpy_vec(yL, m_e[is], vL, NCD);
    }

    int is = m_Ns-1;
    copy_vec(xL, vL, NCD);
    load_vec(vL, &wp2[Nin4*is], NCD);
    add_aPp5_dirac_vec(vL, real_t(0.5)*m_dm[is]/m_dp[is-1], xL, NC);
    add_aPm5_dirac_vec(vL, real_t(-0.5), yL, NC);
    save_vec(&vp2[Nin4*is], vL, NCD);

  }

  #pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::U_inv(AFIELD& v,
                                                  const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t xL[NCD], yL[NCD], vL[NCD];

    int is = m_Ns-1;
    load_vec(vL, &wp2[Nin4*is], NCD);
    scal_vec(vL, real_t(1.0)/(m_dp[is]+m_g), NCD);
    save_vec(&vp2[Nin4*is], vL, NCD);

    set_aPp5_dirac_vec(yL, real_t(0.5), vL, NC);

    for(int is = m_Ns-2; is >= 0; --is){
      copy_vec(xL, vL, NCD);
      load_vec(vL, &wp2[Nin4*is], NCD);
      add_aPm5_dirac_vec(vL, real_t(0.5)*m_dm[is], xL, NC);
      axpy_vec(vL, -m_f[is], yL, NCD);
      scal_vec(vL, real_t(1.0)/m_dp[is], NCD);
      save_vec(&vp2[Nin4*is], vL, NCD);
    }

  }

  #pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Udag_inv(AFIELD& v,
                                                     const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t xL[NCD], yL[NCD], vL[NCD];

    load_vec(vL, &wp2[0], NCD);
    scal_vec(vL, real_t(1.0)/m_dp[0], NCD);
    save_vec(&vp2[0], vL, NCD);
    set_vec(yL, m_f[0], vL, NCD);

    for(int is = 1; is < m_Ns-1; ++is){
      copy_vec(xL, vL, NCD);
      load_vec(vL, &wp2[Nin4*is], NCD);
      add_aPm5_dirac_vec(vL, real_t(0.5)*m_dm[is-1], xL, NC);
      scal_vec(vL, real_t(1.0)/m_dp[is], NCD);
      save_vec(&vp2[Nin4*is], vL, NCD);
      axpy_vec(yL, m_f[is], vL, NCD);
    }

    int is = m_Ns - 1;
    copy_vec(xL, vL, NCD);
    load_vec(vL, &wp2[Nin4*is], NCD);
    add_aPm5_dirac_vec(vL, real_t(0.5)*m_dm[is-1], xL, NC);
    add_aPp5_dirac_vec(vL, real_t(-0.5), yL, NC);
    scal_vec(vL, real_t(1.0)/(m_dp[is] + m_g), NCD);
    save_vec(&vp2[Nin4*is], vL, NCD);

  }

  #pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_General_5din<AFIELD>::Ldag_inv(AFIELD& v,
                                                     const AFIELD& w)
{
  int Nin4 = VLEN * NC * ND;
  int Nin5 = Nin4 * m_Ns;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD*>(&w)->ptr(0);

  int ith, nth, site0, site1;
  set_threadtask_afopr(ith, nth, site0, site1, m_Nstv);

  for(int site = site0; site < site1; ++site){

    real_t *vp2 = &vp[Nin5 * site];
    real_t *wp2 = &wp[Nin5 * site];

    Vsimd_t xL[NCD], yL[NCD], vL[NCD];

    int is = m_Ns-1;
    load_vec(vL, &wp2[Nin4*is], NCD);
    save_vec(&vp2[Nin4*is], vL, NCD);
    set_aPm5_dirac_vec(yL, real_t(0.5), vL, NC);

    for(int is = m_Ns-2; is >= 0; --is){
      copy_vec(xL, vL, NCD);
      load_vec(vL, &wp2[Nin4*is], NCD);
      add_aPp5_dirac_vec(vL, real_t(0.5)*m_dm[is+1]/m_dp[is], xL, NC);
      axpy_vec(vL, -m_e[is], yL, NCD);
      save_vec(&vp2[Nin4*is], vL, NCD);
    }

  }

  #pragma omp barrier

}

//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_General_5din<AFIELD>::flop_count(std::string mode)
{
  int Lvol = CommonParameters::Lvol();
  double vsite = static_cast<double>(Lvol);
  double vNs   = static_cast<double>(m_Ns);

  // double flop_Wilson = m_foprw->flop_count("D");
  double flop_Wilson = m_foprw->flop_count();

  double axpy1 = static_cast<double>(2 * m_NinF);
  double scal1 = static_cast<double>(1 * m_NinF);

  double flop_DW = vNs * (flop_Wilson + vsite*(6*axpy1 + 2*scal1));
  // In Ddag case, flop_Wilson + 7 axpy which equals flop_DW. 

  double flop_LU_inv = 2.0 * vsite * 
            ((3.0*axpy1 + scal1)*(vNs-1.0) + axpy1 + 2.0*scal1);

  double flop = 0.0;
  if(mode == "Prec"){
    flop = flop_LU_inv;
  } else if ((mode == "D") || (mode == "Ddag")){
    flop = flop_DW;
  } else if (mode == "DdagD"){
    flop = 2.0 * flop_DW;
  } else if ((mode == "D_prec") || (mode == "Ddag_prec")){
    flop = flop_LU_inv + flop_DW;
  } else if (mode == "DdagD_prec"){
    flop = 2.0 * (flop_LU_inv + flop_DW);
  }else{
    vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  return flop;

}

//============================================================END=====
