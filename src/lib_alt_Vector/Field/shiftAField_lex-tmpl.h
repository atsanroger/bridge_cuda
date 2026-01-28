/*!
      @file    shiftAField_lex-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

template<typename AFIELD>
const std::string ShiftAField_lex<AFIELD>::class_name =
                                         "ShiftAField_lex<AFIELD>";
//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::init(int Nin)
{
  int Ndim = CommonParameters::Ndim();
  std::vector<int> bc(Ndim);
  for(int mu = 0; mu < Ndim; ++mu) bc[mu] = 1;

  init(Nin, bc);

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::init(int Nin, std::vector<int>& bc)
{
  m_vl = CommonParameters::Vlevel();

  int req_comm = 1;  // set 1 if communication forced any time
  //int req_comm = 0;  // set 0 if communication only when necessary

  vout.general(m_vl, "%s: being constructed.\n", class_name.c_str());

  m_Nin = Nin;
  vout.general(m_vl, "  Nin = %d\n", m_Nin);

  m_Nx = CommonParameters::Nx();
  m_Ny = CommonParameters::Ny();
  m_Nz = CommonParameters::Nz();
  m_Nt = CommonParameters::Nt();
  m_Nvol = m_Nx * m_Ny * m_Nz * m_Nt;

  m_Ndim = CommonParameters::Ndim();

  m_Nsize[0] = m_Nx;
  m_Nsize[1] = m_Ny;
  m_Nsize[2] = m_Nz;
  m_Nsize[3] = m_Nt;

  if(bc.size() != m_Ndim){
    vout.crucial(m_vl, "%s: incorrect size of boundary condition\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_boundary.resize(m_Ndim);

  for(int mu = 0; mu < m_Ndim; ++mu){
    m_boundary[mu] = bc[mu];
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

  do_comm_any = 0;
  for(int mu = 0; mu < m_Ndim; ++mu){
    do_comm[mu] = 1;
    if(req_comm == 0 && Communicator::npe(mu) == 1) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general(m_vl, "  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_bc[mu]  = 1;
    if(do_comm[mu] > 0){ // do communication
      if(Communicator::ipe(mu) == 0) m_bc[mu] = m_boundary[mu];
      m_bc2[mu] = 0;
    }else{  // no communication
      m_bc[mu]  = 0;    // for boundar part (dummy)
      m_bc2[mu] = m_boundary[mu];  // for bulk part
    }
  }



  m_Nbdsize.resize(m_Ndim);
  m_Nbdsize[0] = m_Nin * m_Ny * m_Nz * m_Nt;
  m_Nbdsize[1] = m_Nin * m_Nx * m_Nz * m_Nt;
  m_Nbdsize[2] = m_Nin * m_Nx * m_Ny * m_Nt;
  m_Nbdsize[3] = m_Nin * m_Nx * m_Ny * m_Nz;

  setup_channels();

  vout.general(m_vl, "%s: construction finished.\n", class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::setup_channels()
{
  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for(int mu = 0; mu < m_Ndim; ++mu){

    int Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu,  1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu,  1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void* buf_up = (void*)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu,  1, buf_up);
    void* buf_dn = (void*)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
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
void ShiftAField_lex<AFIELD>::backward(AFIELD& v, const AFIELD& w,
                                                         const int mu)
{
  int Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(v.check_size(m_Nin, m_Nvol, Nex));

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for(int ex = 0; ex < Nex; ++ex){

    real_t* vp = v.ptr(index.idx(0, m_Nin, 0, ex));
    real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, m_Nin, 0, ex));

    if(mu == 0) {
      up_x(vp, wp);
    }else if(mu == 1) {
      up_y(vp, wp);
    }else if(mu == 2) {
      up_z(vp, wp);
    }else if(mu == 3) {
      up_t(vp, wp);
    }else{
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::backward(AFIELD& v, const int ex1,
                                 const AFIELD& w, const int ex2,
                                                  const int mu)
{
  int Nex = v.nex();
  assert(v.check_size(m_Nin, m_Nvol, Nex));
  assert(ex1 < Nex);
  Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(ex2 < Nex);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  real_t* vp = v.ptr(index.idx(0, m_Nin, 0, ex1));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, m_Nin, 0, ex2));

  if(mu == 0) {
    up_x(vp, wp);
  }else if(mu == 1) {
    up_y(vp, wp);
  }else if(mu == 2) {
    up_z(vp, wp);
  }else if(mu == 3) {
    up_t(vp, wp);
  }else{
    vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::forward(AFIELD& v, const AFIELD& w,
                                                         const int mu)
{
  int Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(v.check_size(m_Nin, m_Nvol, Nex));

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for(int ex = 0; ex < Nex; ++ex){

    real_t* vp = v.ptr(index.idx(0, m_Nin, 0, ex));
    real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, m_Nin, 0, ex));

    if (mu == 0) {
      dn_x(vp, wp);
    } else if (mu == 1) {
      dn_y(vp, wp);
    } else if (mu == 2) {
      dn_z(vp, wp);
    } else if (mu == 3) {
      dn_t(vp, wp);
    } else {
      vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::forward(AFIELD& v, const int ex1,
                                const AFIELD& w, const int ex2,
                                                 const int mu)
{
  int Nex = v.nex();
  assert(v.check_size(m_Nin, m_Nvol, Nex));
  assert(ex1 < Nex);
  Nex = w.nex();
  assert(w.check_size(m_Nin, m_Nvol, Nex));
  assert(ex2 < Nex);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  real_t* vp = v.ptr(index.idx(0, m_Nin, 0, ex1));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, m_Nin, 0, ex2));

  if (mu == 0) {
    dn_x(vp, wp);
  } else if (mu == 1) {
    dn_y(vp, wp);
  } else if (mu == 2) {
    dn_z(vp, wp);
  } else if (mu == 3) {
    dn_t(vp, wp);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong parameter\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_x(real_t *v2, real_t *v1)
{
  int idir = 0;

  real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_up[idir].ptr();

  if(do_comm[idir] > 0){
    BridgeVec::shift_lex_xp1(buf1, v1, m_Nin, m_Nsize, m_bc);
    //copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_dn[idir].start();
    chrecv_up[idir].start();

    chsend_dn[idir].wait();
    chrecv_up[idir].wait();
  }

  BridgeVec::shift_lex_xpb(v2, v1, m_Nin, m_Nsize, m_bc2);

  if(do_comm[idir] > 0){
    //  copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeVec::shift_lex_xp2(v2, buf2, m_Nin, m_Nsize, m_bc);
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_x(real_t *v2, real_t *v1)
{
  int idir = 0;

  real_t *buf1 = (real_t*)chsend_up[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();

  if(do_comm[idir] > 0){

    BridgeVec::shift_lex_xm1(buf1, v1, m_Nin, m_Nsize, m_bc);
    // copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_up[idir].start();
    chrecv_dn[idir].start();

    chsend_up[idir].wait();
    chrecv_dn[idir].wait();

  }

  BridgeVec::shift_lex_xmb(v2, v1, m_Nin, m_Nsize, m_bc2);

  if(do_comm[idir] > 0){
    // copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeVec::shift_lex_xm2(v2, buf2, m_Nin, m_Nsize, m_bc);
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_y(real_t *v2, real_t *v1)
{
  int idir = 1;

  real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_up[idir].ptr();

  if(do_comm[idir] > 0){

    BridgeVec::shift_lex_yp1(buf1, v1, m_Nin, m_Nsize, m_bc);
    // copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_dn[idir].start();
    chrecv_up[idir].start();

    chsend_dn[idir].wait();
    chrecv_up[idir].wait();

  }

  BridgeVec::shift_lex_ypb(v2, v1, m_Nin, m_Nsize, m_bc2);

  if(do_comm[idir] > 0){
    // copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeVec::shift_lex_yp2(v2, buf2, m_Nin, m_Nsize, m_bc);
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_y(real_t *v2, real_t *v1)
{
  int idir = 1;

  real_t *buf1 = (real_t*)chsend_up[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();

  if(do_comm[idir] > 0){

    BridgeVec::shift_lex_ym1(buf1, v1, m_Nin, m_Nsize, m_bc);
    // copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_up[idir].start();
    chrecv_dn[idir].start();

    chsend_up[idir].wait();
    chrecv_dn[idir].wait();

  }

  BridgeVec::shift_lex_ymb(v2, v1, m_Nin, m_Nsize, m_bc2);

  if(do_comm[idir] > 0){
    // copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeVec::shift_lex_ym2(v2, buf2, m_Nin, m_Nsize, m_bc);
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_z(real_t *v2, real_t *v1)
{
  int idir = 2;

  real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_up[idir].ptr();

  if(do_comm[idir] > 0){

    BridgeVec::shift_lex_zp1(buf1, v1, m_Nin, m_Nsize, m_bc);
    // copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_dn[idir].start();
    chrecv_up[idir].start();

    chsend_dn[idir].wait();
    chrecv_up[idir].wait();

  }

  BridgeVec::shift_lex_zpb(v2, v1, m_Nin, m_Nsize, m_bc2);

  if(do_comm[idir] > 0){
    // copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeVec::shift_lex_zp2(v2, buf2, m_Nin, m_Nsize, m_bc);
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_z(real_t *v2, real_t *v1)
{
  int idir = 2;

  real_t *buf1 = (real_t*)chsend_up[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();

  if(do_comm[idir] > 0){

    BridgeVec::shift_lex_zm1(buf1, v1, m_Nin, m_Nsize, m_bc);
    // copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_up[idir].start();
    chrecv_dn[idir].start();

    chsend_up[idir].wait();
    chrecv_dn[idir].wait();

  }

  BridgeVec::shift_lex_zmb(v2, v1, m_Nin, m_Nsize, m_bc2);

  if(do_comm[idir] > 0){
    // copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeVec::shift_lex_zm2(v2, buf2, m_Nin, m_Nsize, m_bc);
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::up_t(real_t *v2, real_t *v1)
{
  int idir = 3;

  real_t *buf1 = (real_t*)chsend_dn[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_up[idir].ptr();

  if(do_comm[idir] > 0){

    BridgeVec::shift_lex_tp1(buf1, v1, m_Nin, m_Nsize, m_bc);
    // copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_dn[idir].start();
    chrecv_up[idir].start();

    chsend_dn[idir].wait();
    chrecv_up[idir].wait();

  }

  BridgeVec::shift_lex_tpb(v2, v1, m_Nin, m_Nsize, m_bc2);

  if(do_comm[idir] > 0){
    // copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeVec::shift_lex_tp2(v2, buf2, m_Nin, m_Nsize, m_bc);
  }

}

//====================================================================
template<typename AFIELD>
void ShiftAField_lex<AFIELD>::dn_t(real_t *v2, real_t *v1)
{
  int idir = 3;

  real_t *buf1 = (real_t*)chsend_up[idir].ptr();
  real_t *buf2 = (real_t*)chrecv_dn[idir].ptr();

  if(do_comm[idir] > 0){

    BridgeVec::shift_lex_tm1(buf1, v1, m_Nin, m_Nsize, m_bc);
    // copy_from_device(buf1, m_Nbdsize[idir]);

    chsend_up[idir].start();
    chrecv_dn[idir].start();

    chsend_up[idir].wait();
    chrecv_dn[idir].wait();

  }

  BridgeVec::shift_lex_tmb(v2, v1, m_Nin, m_Nsize, m_bc2);

  if(do_comm[idir] > 0){
    // copy_to_device(buf2, m_Nbdsize[idir]);
    BridgeVec::shift_lex_tm2(v2, buf2, m_Nin, m_Nsize, m_bc);
  }

}

//============================================================END=====
