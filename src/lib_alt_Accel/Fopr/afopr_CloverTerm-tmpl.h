/*!
      @file    afopr_CloverTerm-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-09-10 17:39:07 #$
      @version $LastChangedRevision: 2537 $
*/

template<typename AFIELD>
const std::string AFopr_CloverTerm<AFIELD>::class_name
                                         = "AFopr_CloverTerm<AFIELD>";
//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

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

  m_Nsize[0] = m_Nx;
  m_Nsize[1] = m_Ny;
  m_Nsize[2] = m_Nz;
  m_Nsize[3] = m_Nt;

  //  check_Nc(m_Nc);
  //  check_setup();


  set_parameters(params);

  m_staple = new AStaple_lex<AFIELD>;

  // gauge configuration.
  m_U.reset(m_Ndf, m_Nst, m_Ndim);

  // clover term.
  m_T.reset(m_Ndf, m_Nst, m_Ndm2);
  m_Tinv.reset(m_Ndf, m_Nst, m_Ndm2);

  // working vectors.
  int NinF = 2 * m_Nc * m_Nd;
  m_v1.reset(NinF, m_Nst, 1);
  m_v2.reset(NinF, m_Nst, 1);

  // working gauge field.
  m_ut1.reset(m_Ndf, m_Nst, 1), m_ut2.reset(m_Ndf, m_Nst, 1);

  m_U2.reset(m_Ndf, m_Nst, 4);
  m_F2.reset(m_Ndf, m_Nst, 1);
 
  // setup solver.
  double ecrit = 1.0e-30;
  if(sizeof(real_t) == 4) ecrit = 1.0e-16;

  Parameters params_solver;
  params_solver.set_string("solver_type", "CG");
  params_solver.set_int("maximum_number_of_iteration", 100);
  params_solver.set_int("maximum_number_of_restart", 40);
  params_solver.set_double("convergence_criterion_squared", ecrit);
  //- NB. set VerboseLevel to CRUCIAL to suppress frequent messages.          
  //params_solver.set_string("verbose_level", "Crucial");
  //params_solver.set_string("verbose_level", "Detailed");
  params_solver.set_string("verbose_level", vlevel);

  m_solver = new ASolver_CG<AFIELD>(this);
  m_solver->set_parameters(params_solver);

  vout.decrease_indent();
  vout.detailed(m_vl, "%s: initalization finished.\n",
                      class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_staple;
  delete m_solver;

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

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

  //- setting gamma matrix representation
  std::string repr;
  err = params.fetch_string("gamma_matrix_type", repr);
  if(err){
    vout.general(m_vl, "  gamma_matrix_type is not given - set to Dirac\n");
    m_repr = DIRAC;
  }else if(repr == "Dirac"){
    m_repr = DIRAC;
  }else if(repr == "Chiral"){
    m_repr = CHIRAL;
  }else{
    vout.crucial(m_vl, "Error in %s: irrelevant gamma_matrix_type: %s\n",
		 class_name.c_str(), repr.c_str());
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
  assert(bc.size() == m_Ndim);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_CKs = CKs;
    m_cSW = cSW;
    m_boundary.resize(m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }
  }

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  if(m_repr == DIRAC){
    vout.general(m_vl, "  gamma-matrix type = Dirac\n");
  }else{
    vout.general(m_vl, "  gamma-matrix type = Chiral\n");
  }
  vout.general(m_vl, "  kappa = %8.4f\n", m_CKs);
  vout.general(m_vl, "  cSW   = %8.4f\n", m_cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_config(Field* u)
{
  int nth = ThreadManager::get_num_threads();
  if (nth > 1) {
    set_config_impl(u);
  } else {
    set_config_omp(u);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_config_omp(Field* u)
{
#pragma omp parallel
 {
  set_config_impl(u);
 }
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_config_impl(Field* u)
{
#pragma omp barrier

  vout.general(m_vl, "%s: set_config started.\n", class_name.c_str());

  Timer timer;
  double elapsed_time;
  timer.start();

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_conf = u;

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;

  convert_gauge(index_lex, m_U, *u);

  timer.stop();
  elapsed_time = timer.elapsed_sec();
  vout.general(m_vl, "  convert: %11.6f [sec]\n", elapsed_time);

  timer.reset();
  timer.start();

  set_csw(*u);

  timer.stop();
  elapsed_time = timer.elapsed_sec();
  vout.general(m_vl, "  set_csw: %11.6f [sec]\n", elapsed_time);

  timer.reset();
  timer.start();

  solve_csw_inv();

  timer.stop();
  elapsed_time = timer.elapsed_sec();
  vout.general(m_vl, "  set_csw_inv: %11.6f [sec]\n", elapsed_time);

  vout.general(m_vl, "%s: set_config finished.\n", class_name.c_str());

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_csw(Field& U)
{
  if(m_repr == DIRAC){
    set_csw_dirac(U);
  }else if(m_repr == CHIRAL){
    set_csw_chiral(U);
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

  if(m_repr == DIRAC){
    solve_csw_inv_dirac();
  }else{
    solve_csw_inv_chiral();
  }

  vout.paranoiac(m_vl, "  %s: solving inverse of clover term finished.\n",
                 class_name.c_str());

  // check
  /*
  vout.general(m_vl, "  %s: check of clover term inverse\n",
                 class_name.c_str());

  int NinF = 2 * m_Nc * m_Nd;
  AFIELD w(NinF, m_Nst, 1), w2(NinF, m_Nst, 1);
  AFIELD w3(NinF, m_Nst, 1);

  AIndex_lex<real_t> index;

  for(int id = 0; id < m_Nd; ++id){
    for(int ic = 0; ic < m_Nc; ++ic){

      w.set_host(0.0);
      for(int site = 0; site < m_Nst; ++site){
        w.set_host(index.idx_SPr(ic, id, site, 0), real_t(1.0));
      }
      w.update_device();

      mult_csw(w3, w);
      mult_csw_inv(w2, w3);
      axpy(w2, real_t(-1.0), w);
      real_t ww = w2.norm2();

      vout.general(m_vl, "  ic = %d id = %d  diff2 = %f\n",
		   ic, id, ww);

    }
  }
  */

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::solve_csw_inv_dirac()
{
#pragma omp barrier

  set_mode("D");

  int Nconv;
  real_t diff;
  AIndex_lex<real_t,AFIELD::IMPL> index;

  m_Tinv.set_host(0.0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  int Nd2  = m_Nd/2;
  for(int id = 0; id < Nd2; ++id){
    for(int ic = 0; ic < m_Nc; ++ic){

      m_v1.set_host(0.0);
      for(int site = is; site < ns; ++site){
        m_v1.set_host(index.idx_SPr(ic, id, site, 0), real_t(1.0));
      }
#pragma omp barrier

      update_device(m_v1);
      copy(m_v2, m_v1);

      m_solver->solve(m_v2, m_v1, Nconv, diff);
      vout.paranoiac(m_vl, "  ic = %d  id = %d  Nconv = %d  diff = %12.4e\n",
                     ic, id, Nconv, diff);

      update_host(m_v2);

      for(int site = is; site < ns; ++site){
        for(int ic2 = 0; ic2 < m_Nc; ++ic2){
          for(int id2 = 0; id2 < m_Nd; ++id2){
            real_t re = m_v2.cmp_host(index.idx_SPr(ic2, id2, site, 0));
            real_t im = m_v2.cmp_host(index.idx_SPi(ic2, id2, site, 0));
            int iT = id2 + m_Nd * id;
            m_Tinv.set_host(index.idx_Gr(ic2, ic, site, iT),  re);
            m_Tinv.set_host(index.idx_Gi(ic2, ic, site, iT), -im);
          }
        }
      }
#pragma omp barrier

    }
  }

  update_device(m_Tinv);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::solve_csw_inv_chiral()
{
#pragma omp barrier

  set_mode("D");

  int Nconv;
  real_t diff;
  AIndex_lex<real_t,AFIELD::IMPL> index;

  m_Tinv.set_host(0.0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nst);

  int Nd2  = m_Nd/2;
  for(int id = 0; id < Nd2; ++id){
    for(int ic = 0; ic < m_Nc; ++ic){

      m_v1.set_host(0.0);
      for(int site = is; site < ns; ++site){
        m_v1.set_host(index.idx_SPr(ic, id, site, 0), real_t(1.0));
        m_v1.set_host(index.idx_SPr(ic, id+Nd2, site, 0), real_t(1.0));
      }
#pragma omp barrier

      update_device(m_v1);
      copy(m_v2, m_v1);

      m_solver->solve(m_v2, m_v1, Nconv, diff);
      vout.paranoiac(m_vl, "  ic = %d  id = %d  Nconv = %d  diff = %12.4e\n",
                     ic, id, Nconv, diff);

      update_host(m_v2);

      for(int site = is; site < ns; ++site){

        for(int ic2 = 0; ic2 < m_Nc; ++ic2){
          for(int id2 = 0; id2 < Nd2; ++id2){
            real_t re = m_v2.cmp_host(index.idx_SPr(ic2, id2, site, 0));
            real_t im = m_v2.cmp_host(index.idx_SPi(ic2, id2, site, 0));
            int iT = id2 + Nd2 * id;
            m_Tinv.set_host(index.idx_Gr(ic2, ic, site, iT),  re);
            m_Tinv.set_host(index.idx_Gi(ic2, ic, site, iT), -im);
          }
        }

        for(int ic2 = 0; ic2 < m_Nc; ++ic2){
          for(int id2 = 0; id2 < Nd2; ++id2){
            int jd2 = id2 + Nd2;
            real_t re = m_v2.cmp_host(index.idx_SPr(ic2, jd2, site, 0));
            real_t im = m_v2.cmp_host(index.idx_SPi(ic2, jd2, site, 0));
            int iT = id2 + Nd2 * id + ND;
            m_Tinv.set_host(index.idx_Gr(ic2, ic, site, iT),  re);
            m_Tinv.set_host(index.idx_Gi(ic2, ic, site, iT), -im);
          }
        }

      }
#pragma omp barrier

    }
  }

  update_device(m_Tinv);

#pragma omp barrier

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

  // the following is temporary prescription.
  // to be performed in set_config.   [H.Matsufuru 16 May 2021]
  Timer timer;

  timer.start();

  solve_csw_inv();

  timer.stop();
  double elapsed_time = timer.elapsed_sec();
  vout.detailed(m_vl, "%s: set_csw_inv: %11.6f [sec]\n",
               class_name.c_str(), elapsed_time);


  //#pragma omp parallel
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

#pragma omp barrier
  //  AFIELD U2(m_Ndf, m_Nst, 4), F2(m_Ndf, m_Nst, 1);

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;
  //  convert_gauge(index_lex, U2, U);

  //#pragma omp parallel
 {
  convert_gauge(index_lex, m_U2, U);

  m_T.set(0.0);
  #pragma omp barrier

  //- sigma23
  set_fieldstrength(m_F2, m_U2, 1, 2);
  xI(m_F2);

#pragma omp barrier
  axpy(m_T, 1,  real_t(1.0), m_F2, 0);
  axpy(m_T, 4,  real_t(1.0), m_F2, 0);
#pragma omp barrier

  //- sigma31
  set_fieldstrength(m_F2, m_U2, 2, 0);
#pragma omp barrier
  axpy(m_T, 1, real_t( 1.0), m_F2, 0);
  axpy(m_T, 4, real_t(-1.0), m_F2, 0);
#pragma omp barrier

  //- sigma12
  set_fieldstrength(m_F2, m_U2, 0, 1);
  xI(m_F2);
#pragma omp barrier
  axpy(m_T, 0, real_t( 1.0), m_F2, 0);
  axpy(m_T, 5, real_t(-1.0), m_F2, 0);
#pragma omp barrier

  //- sigma41
  set_fieldstrength(m_F2, m_U2, 3, 0);
  xI(m_F2);
#pragma omp barrier
  axpy(m_T, 3, real_t(-1.0), m_F2, 0);
  axpy(m_T, 6, real_t(-1.0), m_F2, 0);
#pragma omp barrier

  //- sigma42
  set_fieldstrength(m_F2, m_U2, 3, 1);
#pragma omp barrier
  axpy(m_T, 3, real_t(-1.0), m_F2, 0);
  axpy(m_T, 6, real_t( 1.0), m_F2, 0);
#pragma omp barrier

  //- sigma43
  set_fieldstrength(m_F2, m_U2, 3, 2);
  xI(m_F2);
#pragma omp barrier
  axpy(m_T, 2, real_t(-1.0), m_F2, 0);
  axpy(m_T, 7, real_t( 1.0), m_F2, 0);
#pragma omp barrier

  scal(m_T, -m_CKs * m_cSW);

#pragma omp barrier

  add_unit(m_T, 0, real_t(1.0));
  add_unit(m_T, 5, real_t(1.0));

#pragma omp barrier
 }

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::set_csw_chiral(Field& U)
{ //  this is for Dirac representation.

  // The clover term in the chiral representation is as spin-space
  // matrix
  //   [ P+Q   0  ]
  //   [  0   P-Q ],
  // where P and Q are 2x2 block matrices as
  //   P =  [          iF(1,2)   F(3,1) + iF(2,3) ]
  //        [-F(3,1) + iF(2,3)          - iF(1,2) ]
  // and
  //   Q =  [        - iF(4,3)  -F(4,2) - iF(4,1) ]
  //        [ F(4,2) - iF(4,1)            iF(4,3) ]
  // up to the coefficient.
  // in the following what defined is
  //  [ T(0) T(1) ] = P + Q   [ T(4) T(5) ] = P - Q
  //  [ T(2) T(3) ]           [ T(6) T(7) ]
  // T = 1 - kappa c_SW sigma F / 2

  //  AFIELD U2(m_Ndf, m_Nst, 4), F2(m_Ndf, m_Nst, 1);

#pragma omp barrier

  AIndex_lex<real_t,AFIELD::IMPL> index_lex;
  //  convert_gauge(index_lex, U2, U);

  //#pragma omp parallel
 {
  convert_gauge(index_lex, m_U2, U);

  m_T.set(0.0);
  #pragma omp barrier

  //- sigma23
  set_fieldstrength(m_F2, m_U2, 1, 2);  // F_23
  xI(m_F2);
#pragma omp barrier
  axpy(m_T, 1,  real_t(1.0), m_F2, 0);
  axpy(m_T, 2,  real_t(1.0), m_F2, 0);
  axpy(m_T, 5,  real_t(1.0), m_F2, 0);
  axpy(m_T, 6,  real_t(1.0), m_F2, 0);
#pragma omp barrier

  //- sigma31
  set_fieldstrength(m_F2, m_U2, 2, 0);  // F_31
#pragma omp barrier
  axpy(m_T, 1, real_t( 1.0), m_F2, 0);
  axpy(m_T, 2, real_t(-1.0), m_F2, 0);
  axpy(m_T, 5, real_t( 1.0), m_F2, 0);
  axpy(m_T, 6, real_t(-1.0), m_F2, 0);
#pragma omp barrier

  //- sigma12
  set_fieldstrength(m_F2, m_U2, 0, 1);  // F_12
  xI(m_F2);
#pragma omp barrier
  axpy(m_T, 0, real_t( 1.0), m_F2, 0);
  axpy(m_T, 3, real_t(-1.0), m_F2, 0);
  axpy(m_T, 4, real_t( 1.0), m_F2, 0);
  axpy(m_T, 7, real_t(-1.0), m_F2, 0);
#pragma omp barrier

  //- sigma41
  set_fieldstrength(m_F2, m_U2, 3, 0);  // F_41
  xI(m_F2);
#pragma omp barrier
  axpy(m_T, 1, real_t(-1.0), m_F2, 0);
  axpy(m_T, 2, real_t(-1.0), m_F2, 0);
  axpy(m_T, 5, real_t( 1.0), m_F2, 0);
  axpy(m_T, 6, real_t( 1.0), m_F2, 0);
#pragma omp barrier

  //- sigma42
  set_fieldstrength(m_F2, m_U2, 3, 1);  // F_42
#pragma omp barrier
  axpy(m_T, 1, real_t(-1.0), m_F2, 0);
  axpy(m_T, 2, real_t( 1.0), m_F2, 0);
  axpy(m_T, 5, real_t( 1.0), m_F2, 0);
  axpy(m_T, 6, real_t(-1.0), m_F2, 0);
#pragma omp barrier

  //- sigma43
  set_fieldstrength(m_F2, m_U2, 3, 2);  // F_43
  xI(m_F2);
#pragma omp barrier
  axpy(m_T, 0, real_t(-1.0), m_F2, 0);
  axpy(m_T, 3, real_t( 1.0), m_F2, 0);
  axpy(m_T, 4, real_t( 1.0), m_F2, 0);
  axpy(m_T, 7, real_t(-1.0), m_F2, 0);
#pragma omp barrier

  scal(m_T, -m_CKs * m_cSW);

#pragma omp barrier

  add_unit(m_T, 0, real_t(1.0));
  add_unit(m_T, 3, real_t(1.0));
  add_unit(m_T, 4, real_t(1.0));
  add_unit(m_T, 7, real_t(1.0));

#pragma omp barrier
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
void AFopr_CloverTerm<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
std::string AFopr_CloverTerm<AFIELD>::get_mode() const
{
  return m_mode;
}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult(AFIELD &v, const AFIELD &w)
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
void AFopr_CloverTerm<AFIELD>::mult_dag(AFIELD &v, const AFIELD &w)
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
#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    if(m_repr == DIRAC){
      BridgeACC::mult_wilson_gm5_dirac(v, w, m_Nsize, NC);
    }else{
      BridgeACC::mult_wilson_gm5_chiral(v, w, m_Nsize, NC);
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_csw(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    real_t *v2 = v.ptr(0);
    real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
    real_t *u  = m_T.ptr(0);

    if(m_repr == DIRAC){
      BridgeACC::mult_csw_dirac(v2, u, v1, m_Nsize, 0, NC);
    }else{
      BridgeACC::mult_csw_chiral(v2, u, v1, m_Nsize, 0, NC);
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::mult_csw_inv(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *u  = m_Tinv.ptr(0);

#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    if(m_repr == DIRAC){
      BridgeACC::mult_csw_dirac(v2, u, v1, m_Nsize, 0, NC);
    }else{
      BridgeACC::mult_csw_chiral(v2, u, v1, m_Nsize, 0, NC);
    }
  }

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_CloverTerm<AFIELD>::multadd_csw(AFIELD& v, const AFIELD& w)
{
  real_t *v2 = v.ptr(0);
  real_t *v1 = const_cast<AFIELD*>(&w)->ptr(0);
  real_t *u  = m_T.ptr(0);

#pragma omp barrier

  int ith, nth;
  set_thread(ith, nth);

  if(ith == 0){
    if(m_repr == DIRAC){
      BridgeACC::mult_csw_dirac(v2, u, v1, m_Nsize, 1, NC);
    }else{
      BridgeACC::mult_csw_chiral(v2, u, v1, m_Nsize, 1, NC);
    }
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

  if (m_repr == DIRAC) {
    flop_site = static_cast<double>(
      m_Nc * m_Nd * (4 + 6 * (4 * m_Nc + 2) + 2 * (4 * m_Nc + 1))
      + 8 * m_Nc * m_Nc * m_Nd * m_Nd); // <- clover term
  } else if (m_repr == CHIRAL) {
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
