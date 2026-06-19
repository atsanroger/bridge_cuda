/*!
      @file    afopr_CloverTerm-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/


template<typename AFIELD>
const std::string AFopr_CloverTerm<AFIELD>::class_name
                                         = "AFopr_CloverTerm<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::init()
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  m_repr = "Dirac";  // now only the Dirac repr is available.

  // switches
  int req_comm = 1;  // set 1 if communication forced any time
  // int req_comm = 0;  // set 1 if communication forced any time

  m_vl = AFopr<AFIELD>::m_vl;

  m_Nc = CommonParameters::Nc();
  m_Nd = CommonParameters::Nd();
  m_Ndim = CommonParameters::Ndim();
  m_Nvc  = m_Nc * 2;
  m_Ndf  = 2 * m_Nc * m_Nc;
  m_Ndm2 = m_Nd * m_Nd / 2,

  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nst  = CommonParameters::Nvol();

  check_Nc(m_Nc);
  check_setup();

  m_Nxv  = m_Nx/VLEN2;
  m_Nstv = m_Nst/VLEN2;
  vout.general(m_vl, "  Nxv = %d  Nstv = %d\n", m_Nxv, m_Nstv);

  m_staple = new AStaple_lex<AFIELD>;

  // gauge configuration.
  m_U.reset(m_Ndf, m_Nst, m_Ndim);

  // clover term.
  m_T.reset(m_Ndf, m_Nst, m_Ndm2);
  m_Tinv.reset(m_Ndf, m_Nst, m_Ndm2);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v2.reset(NinF, m_Nst, 1);

  // working gauge field.
  m_ut1.reset(m_Ndf, m_Nst, 1), m_ut2.reset(m_Ndf, m_Nst, 1);

  // setup solver.
  Parameters params_solver;
  params_solver.set_string("solver_type", "CG");
  params_solver.set_int("maximum_number_of_iteration", 100);
  params_solver.set_int("maximum_number_of_restart", 40);
  params_solver.set_double("convergence_criterion_squared", 1.0e-30);
  //- NB. set VerboseLevel to CRUCIAL to suppress frequent messages.          
  params_solver.set_string("verbose_level", "Crucial");

  m_solver = new ASolver_CG<AFIELD>(this);
  m_solver->set_parameters(params_solver);


  vout.detailed(m_vl, "%s: initalization finished.\n",
  		      class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::tidyup()
{
  //  ThreadManager_OpenMP::assert_single_thread(class_name);

  delete m_staple;
  delete m_solver;

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double kappa, cSW;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);

  err += params.fetch_int_vector("boundary_condition", bc);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                                                class_name.c_str());
  exit(EXIT_FAILURE);
 }

  set_parameters(real_t(kappa), real_t(cSW), bc);

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_parameters(const real_t CKs,
                                          const real_t cSW,
                                          const std::vector<int> bc)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  assert(bc.size() == m_Ndim);

  m_CKs = CKs;
  m_cSW = cSW;
  m_boundary.resize(m_Ndim);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  CKs  = %8.4f\n", m_CKs);
  vout.general(m_vl, "  cSW  = %8.4f\n", m_cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_config(Field* u)
{
  ThreadManager_OpenMP::assert_single_thread(class_name);

  Timer timer;
  double elapsed_time;
  timer.start();

  m_conf = u;

  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;

#pragma omp parallel
 {
  convert_gauge(index_lex, m_U, *u);
 }

  for(int mu = 0; mu < m_Ndim; ++mu){
    if(m_boundary[mu] != 1) set_boundary_config(m_U, mu);
  }

  timer.stop();
  elapsed_time = timer.elapsed_sec();
  vout.detailed(m_vl, "%s: convert and set boundary: %11.6f [sec]\n",
               class_name.c_str(), elapsed_time);

  timer.reset();
  timer.start();

  set_csw(*u);

  timer.stop();
  elapsed_time = timer.elapsed_sec();
  vout.detailed(m_vl, "%s: set_csw: %11.6f [sec]\n",
               class_name.c_str(), elapsed_time);

  timer.reset();
  timer.start();

  solve_csw_inv();

  timer.stop();
  elapsed_time = timer.elapsed_sec();
  vout.detailed(m_vl, "%s: set_csw_inv: %11.6f [sec]\n",
               class_name.c_str(), elapsed_time);

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_csw(Field& U)
{
  if(m_repr == "Dirac"){
    set_csw_dirac(U);
  }else if(m_repr == "Chiral"){
    // set_csw_chiral(U);
    vout.crucial(m_vl, "%s: Chiral repr has not been suppoted.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);

  }else{
    vout.crucial(m_vl, "%s: unsupported representation.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::solve_csw_inv()
{
  vout.paranoiac(m_vl, "  %s: solving inverse of clover term started.\n",
               class_name.c_str());

  set_mode("D");

  int NinF = 2 * m_Nc * m_Nd;
  AFIELD w(NinF, m_Nst, 1), w2(NinF, m_Nst, 1);

  int Nconv;
  real_t diff;
  Index_lex_alt<real_t,AFIELD::IMPL> index;

#pragma omp parallel
 {
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  // the following code is for the Dirac repr.
  int Nd2  = m_Nd/2;
  for(int id = 0; id < Nd2; ++id){
    for(int ic = 0; ic < m_Nc; ++ic){

      w.set(0.0);

      for(int site = is; site < ns; ++site){
        w.set(index.idx_SPr(ic, id, site, 0), real_t(1.0));
      }

      copy(w2, w);
      m_solver->solve(w2, w, Nconv, diff);
      vout.paranoiac(m_vl, "  ic = %d  id = %d  Nconv = %d  diff = %12.4e\n",
                     ic, id, Nconv, diff);

      for(int site = is; site < ns; ++site){
        for(int ic2 = 0; ic2 < m_Nc; ++ic2){
          for(int id2 = 0; id2 < m_Nd; ++id2){
            real_t re = w2.cmp(index.idx_SPr(ic2, id2, site, 0));
            real_t im = w2.cmp(index.idx_SPi(ic2, id2, site, 0));
            int iT = id2 + m_Nd * id;
            m_Tinv.set(index.idx_Gr(ic2, ic, site, iT),  re);
            m_Tinv.set(index.idx_Gi(ic2, ic, site, iT), -im);
	  }
	}
      }

    }
  }

 } // omp parallel

  vout.paranoiac(m_vl, "  %s: solving inverse of clover term finished.\n",
                 class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::get_csw(AFIELD& T)
{
  if(T.check_size(m_T) == false){
    vout.crucial(m_vl, "%s: in get_csw, incorrect AFIELD size.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  copy(T, m_T);

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::get_csw_inv(AFIELD& T)
{
  if(T.check_size(m_T) == false){
    vout.crucial(m_vl, "%s: in get_csw_inv, incorrect AFIELD size.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp parallel
 {
  copy(T, m_Tinv);
 }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::get_csw_inv(AFIELD& T, const int j)
{
  int Ndf = m_T.nin();
  int Nst = m_T.nvol();
  int Nex = m_T.nex();
  if(T.check_size(Ndf, Nst, 1) == false){
    vout.crucial(m_vl, "%s: in get_csw_inv, incorrect AFIELD size.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if(j >= Nex){
    vout.crucial(m_vl, "%s: in get_csw_inv, index is too large.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

#pragma omp parallel
 {
  copy(T, 0, m_Tinv, j);
 }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_csw_dirac(Field& U)
{ //  this is for Dirac representation.

  // The clover term in the Dirac representation is as spin-space
  // matrix
  //   [ P Q ]
  //   [ Q P ],
  // where P and Q are 2x2 block matrices as
  //   P =  [          iF(1,2)   F(3,1) + iF(2,3) ]
  //        [-F(3,1) + iF(2,3)          - iF(1,2) ]
  // and
  //   Q =  [        - iF(4,3)  -F(4,2) - iF(4,1) ]
  //        [ F(4,2) - iF(4,1)            iF(4,3) ]
  // up to the coefficient.
  // in the following what defined is
  // [ P Q ] = [ T(0) T(1)  T(2) T(3) ]
  //           [ T(4) T(5)  T(6) T(7) ].

  //! T = 1 - kappa c_SW sigma F / 2

  AFIELD U2(m_Ndf, m_Nst, 4), F2(m_Ndf, m_Nst, 1);

  Index_lex_alt<real_t,AFIELD::IMPL> index_lex;
  //  convert_gauge(index_lex, U2, U);

#pragma omp parallel
 {
  convert_gauge(index_lex, U2, U);

  m_T.set(0.0);
  #pragma omp barrier

  //- sigma23
  set_fieldstrength(F2, U2, 1, 2);
  xI(F2);
#pragma omp barrier
  axpy(m_T, 1,  real_t(1.0), F2, 0);
  axpy(m_T, 4,  real_t(1.0), F2, 0);
#pragma omp barrier

  //- sigma31
  set_fieldstrength(F2, U2, 2, 0);
#pragma omp barrier
  axpy(m_T, 1, real_t( 1.0), F2, 0);
  axpy(m_T, 4, real_t(-1.0), F2, 0);
#pragma omp barrier

  //- sigma12
  set_fieldstrength(F2, U2, 0, 1);
  xI(F2);
#pragma omp barrier
  axpy(m_T, 0, real_t( 1.0), F2, 0);
  axpy(m_T, 5, real_t(-1.0), F2, 0);
#pragma omp barrier

  //- sigma41
  set_fieldstrength(F2, U2, 3, 0);
  xI(F2);
#pragma omp barrier
  axpy(m_T, 3, real_t(-1.0), F2, 0);
  axpy(m_T, 6, real_t(-1.0), F2, 0);
#pragma omp barrier

  //- sigma42
  set_fieldstrength(F2, U2, 3, 1);
#pragma omp barrier
  axpy(m_T, 3, real_t(-1.0), F2, 0);
  axpy(m_T, 6, real_t( 1.0), F2, 0);
#pragma omp barrier

  //- sigma43
  set_fieldstrength(F2, U2, 3, 2);
  xI(F2);
#pragma omp barrier
  axpy(m_T, 2, real_t(-1.0), F2, 0);
  axpy(m_T, 7, real_t( 1.0), F2, 0);
#pragma omp barrier

  scal(m_T, -m_CKs * m_cSW);

#pragma omp barrier

  add_unit(m_T, 0, real_t(1.0));
  add_unit(m_T, 5, real_t(1.0));

 }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_fieldstrength(AFIELD& Fst, AFIELD& U,
                                               int mu, int nu)
{
  m_staple->upper_th(m_ut2, U, mu, nu);
  mult_Gnd(Fst, 0, U, mu, m_ut2, 0);
  mult_Gdn(m_ut1, 0, m_ut2, 0, U, mu);

  m_staple->lower_th(m_ut2, U, mu, nu);
  multadd_Gnd(Fst, 0, U, mu, m_ut2, 0, real_t(-1.0));
  multadd_Gdn(m_ut1, 0, m_ut2, 0, U, mu, real_t(-1.0));

  m_staple->shift_forward(m_ut2, 0, m_ut1, 0, mu);

  axpy(Fst, real_t(1.0), m_ut2);
#pragma omp barrier
  ah_G(Fst, 0);
#pragma omp barrier
  scal(Fst, real_t(0.25));
#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_boundary_config(AFIELD& U,
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

#pragma omp parallel
 {
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);
 
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

}

//====================================================================
template<typename AFIELD>
std::string AFopr_CloverTerm<AFIELD>::get_mode() const
{
  return m_mode;
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult(AFIELD &v,
                                  const AFIELD &w)
{
  if(m_mode == "D"){
    D(v, w);
  }else if(m_mode == "Dinv"){
    mult_csw_inv(v, w);
  }else if(m_mode == "H"){
    return H(v, w);
  }else{
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_dag(AFIELD &v,
                                  const AFIELD &w)
{
  if(m_mode == "D"){
    D(v, w);
  }else if(m_mode == "Dinv"){
    mult_csw_inv(v, w);
  }else if(m_mode == "H"){
    H(v, w);
  }else{
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult(AFIELD &v, const AFIELD &w,
                                const std::string mode)
{
  if(mode == "D"){
    D(v, w);
  }else if(mode == "H"){
    H(v, w);
  }else{
    vout.crucial(m_vl, "%s: illegal mode is given to mult with mode\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm5(AFIELD &v, const AFIELD &w)
{
  real_t* vp = v.ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_gm5(vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::D(AFIELD &v, const AFIELD &w)
{
  mult_csw(v, w);
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::H(AFIELD &v, const AFIELD &w)
{
  mult_csw(m_v2, w);
  mult_gm5(v, m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_gm5(real_t *v, real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  for(int site = is; site < ns; ++site){
    for(int ic = 0; ic < m_Nc; ++ic){
      v[index.idx_SPr(ic,0,site,0)] = w[index.idx_SPr(ic,2,site,0)];
      v[index.idx_SPi(ic,0,site,0)] = w[index.idx_SPi(ic,2,site,0)];
      v[index.idx_SPr(ic,1,site,0)] = w[index.idx_SPr(ic,3,site,0)];
      v[index.idx_SPi(ic,1,site,0)] = w[index.idx_SPi(ic,3,site,0)];
      v[index.idx_SPr(ic,2,site,0)] = w[index.idx_SPr(ic,0,site,0)];
      v[index.idx_SPi(ic,2,site,0)] = w[index.idx_SPi(ic,0,site,0)];
      v[index.idx_SPr(ic,3,site,0)] = w[index.idx_SPr(ic,1,site,0)];
      v[index.idx_SPi(ic,3,site,0)] = w[index.idx_SPi(ic,1,site,0)];
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_csw(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);

  Vsimd_t v2v[NCD];
  real_t zL[VLEN*NCD], uL[VLEN*NDF2];

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

#pragma omp barrier

  for(int site = is; site < ns; ++site){

    clear_vec(v2v, NCD);

    //    scal_vec(v2v, -m_CKs, NCD);

    Vsimd_t v1v[NCD];
    load_vec(v1v, &v1[VLEN*NCD*site], NCD);

    // this is in Dirac representation
    Vsimd_t ut[NDF2], wt1, wt2;
    real_t *u = m_T.ptr(0);
    for(int jd = 0; jd < ND2; ++jd){
      for(int id = 0; id < ND; ++id){
        int ig = VLEN * NDF2 * (site + m_Nstv * (id + ND * jd));
        load_vec(ut, &u[ig], NDF2);
        for(int ic = 0; ic < m_Nc; ++ic){
          int ic2 = NC * ic;
          int id2 = (id + ND2) % ND;
          mult_uv(wt1, &ut[ic2], &v1v[NC*id],  NC);
          mult_uv(wt2, &ut[ic2], &v1v[NC*id2], NC);
          int icd1 = ic + NC * jd;
          int icd2 = ic + NC * (jd + ND2);
          axpy_vec(&v2v[icd1], real_t(1.0), &wt1, 1);
          axpy_vec(&v2v[icd2], real_t(1.0), &wt2, 1);
	}
      }
    }

    save_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_csw_inv(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);

  Vsimd_t v2v[NCD];
  real_t zL[VLEN*NCD], uL[VLEN*NDF2];

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

#pragma omp barrier

  for(int site = is; site < ns; ++site){

    clear_vec(v2v, NCD);

    Vsimd_t v1v[NCD];
    load_vec(v1v, &v1[VLEN*NCD*site], NCD);

    // this is in Dirac representation
    Vsimd_t ut[NDF2], wt1, wt2;
    real_t *u = m_Tinv.ptr(0);
    for(int jd = 0; jd < ND2; ++jd){
      for(int id = 0; id < ND; ++id){
        int ig = VLEN * NDF2 * (site + m_Nstv * (id + ND * jd));
        load_vec(ut, &u[ig], NDF2);
        for(int ic = 0; ic < m_Nc; ++ic){
          int ic2 = NC * ic;
          int id2 = (id + ND2) % ND;
          mult_uv(wt1, &ut[ic2], &v1v[NC*id],  NC);
          mult_uv(wt2, &ut[ic2], &v1v[NC*id2], NC);
          int icd1 = ic + NC * jd;
          int icd2 = ic + NC * (jd + ND2);
          axpy_vec(&v2v[icd1], real_t(1.0), &wt1, 1);
          axpy_vec(&v2v[icd2], real_t(1.0), &wt2, 1);
	}
      }
    }

    save_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::multadd_csw(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);

  Vsimd_t v2v[NCD];
  real_t zL[VLEN*NCD], uL[VLEN*NDF2];

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

#pragma omp barrier

  for(int site = is; site < ns; ++site){

    load_vec(v2v, &v2[VLEN*NCD*site], NCD);

    //    scal_vec(v2v, -m_CKs, NCD);

    Vsimd_t v1v[NCD];
    load_vec(v1v, &v1[VLEN*NCD*site], NCD);

    // this is in Dirac representation
    Vsimd_t ut[NDF2], wt1, wt2;
    real_t *u = m_T.ptr(0);
    for(int jd = 0; jd < ND2; ++jd){
      for(int id = 0; id < ND; ++id){
        int ig = VLEN * NDF2 * (site + m_Nstv * (id + ND * jd));
        load_vec(ut, &u[ig], NDF2);
        for(int ic = 0; ic < m_Nc; ++ic){
          int ic2 = NC * ic;
          int id2 = (id + ND2) % ND;
          mult_uv(wt1, &ut[ic2], &v1v[NC*id],  NC);
          mult_uv(wt2, &ut[ic2], &v1v[NC*id2], NC);
          int icd1 = ic + NC * jd;
          int icd2 = ic + NC * (jd + ND2);
          axpy_vec(&v2v[icd1], real_t(1.0), &wt1, 1);
          axpy_vec(&v2v[icd2], real_t(1.0), &wt2, 1);
	}
      }
    }

    save_vec(&v2[VLEN*NCD*site], v2v, NCD);

  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
double AFopr_CloverTerm<AFIELD>::flop_count()
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

  int    Lvol = CommonParameters::Lvol();
  double flop_site, flop;

  if (m_repr == "Dirac") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1))
      + 8 * m_Nc * m_Nc * m_Nd * m_Nd); // <- clover term
  } else if (m_repr == "Chiral") {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 8 * (4 * m_Nc + 2))
      + 8 * m_Nc * m_Nc * m_Nd * m_Nd);   // <- clover term
  } else {
    //    vout.crucial(m_vl, "%s: input repr is undefined.\n",
    //                 class_name.c_str());
    vout.crucial(m_vl, "%s: input repr is undefined.\n");
    abort();
  }

  flop = flop_site * static_cast<double>(Lvol);
  if ((m_mode == "DdagD") || (m_mode == "DDdag")) flop *= 2.0;

  return flop;
}

//============================================================END=====
