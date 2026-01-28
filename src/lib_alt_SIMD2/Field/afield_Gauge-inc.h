/*!
      @file    afield_Gauge-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef AFIELD_GAUGE_INC_INCLUDED
#define AFIELD_GAUGE_INC_INCLUDED

//#include <cstdlib> 


//====================================================================
template <typename REALTYPE>
void multadd_Gnn(AField<REALTYPE,SIMD2>& u, const int exu,
                 const AField<REALTYPE,SIMD2>& v, const int exv,
                 const AField<REALTYPE,SIMD2>& w, const int exw,
                 const REALTYPE a)
{
  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NDF, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NDF, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, ut2, vt[NC], wt[NDF2];

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&wt[ic2+NC*ic1], &wp[iv2], 1);
      }
    }

    for(int ic2 = 0; ic2 < NC; ++ic2){
      int iv1 = iv + VLEN * NC * ic2;
      load_vec(vt, &vp[iv1], NC);
      for(int ic1 = 0; ic1 < NC; ++ic1){
	// mult_uv(ut, vt, &wt[NC*ic1], NC);
        mult_Vnn(ut, vt, &wt[NC*ic1], NC);
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&ut2, &up[iv2], 1);
        axpy_vec(&ut2, a, &ut,1);
        save_vec(&up[iv2], &ut2, 1);
      }
    }

  }

}

//====================================================================
template <typename REALTYPE>
void mult_Gnn(AField<REALTYPE,SIMD2>& u, const int exu,
              const AField<REALTYPE,SIMD2>& v, const int exv,
              const AField<REALTYPE,SIMD2>& w, const int exw)
{
  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NDF, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NDF, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, vt[NC], wt[NDF2];

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&wt[ic2+NC*ic1], &wp[iv2], 1);
      }
    }

    for(int ic2 = 0; ic2 < NC; ++ic2){
      int iv1 = iv + VLEN * NC * ic2;
      load_vec(vt, &vp[iv1], NC);
      for(int ic1 = 0; ic1 < NC; ++ic1){
	// mult_uv(ut, vt, &wt[NC*ic1], NC);
        mult_Vnn(ut, vt, &wt[NC*ic1], NC);
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        save_vec(&up[iv2], &ut, 1);
      }
    }

  }

}

//====================================================================
template <typename REALTYPE>
void multadd_Gnd(AField<REALTYPE,SIMD2>& u, const int exu,
                 const AField<REALTYPE,SIMD2>& v, const int exv,
                 const AField<REALTYPE,SIMD2>& w, const int exw,
                 const REALTYPE a)
{
  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NDF, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NDF, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, ut2, vt[NC], wt[NDF2];

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&wt[ic1+NC*ic2], &wp[iv2], 1);
      }
    }
    // this can be replaced by load_vec(wt, &wp[iv], NDF2);

    for(int ic2 = 0; ic2 < NC; ++ic2){
      int iv1 = iv + VLEN * NC * ic2;
      load_vec(vt, &vp[iv1], NC);
      for(int ic1 = 0; ic1 < NC; ++ic1){
        mult_Vdn(ut, &wt[NC*ic1], vt, NC); // instead of uvdag(vt, wt)
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&ut2, &up[iv2], 1);
        axpy_vec(&ut2, a, &ut, 1);
        save_vec(&up[iv2], &ut2, 1);
      }
    }

  }

}

//====================================================================
template <typename REALTYPE>
void mult_Gnd(AField<REALTYPE,SIMD2>& u, const int exu,
              const AField<REALTYPE,SIMD2>& v, const int exv,
              const AField<REALTYPE,SIMD2>& w, const int exw)
{
  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NDF, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NDF, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, vt[NC], wt[NDF2];

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&wt[ic1+NC*ic2], &wp[iv2], 1);
      }
    }
    // this can be replaced by load_vec(wt, &wp[iv], NDF2);

    for(int ic2 = 0; ic2 < NC; ++ic2){
      int iv1 = iv + VLEN * NC * ic2;
      load_vec(vt, &vp[iv1], NC);
      for(int ic1 = 0; ic1 < NC; ++ic1){
        mult_Vdn(ut, &wt[NC*ic1], vt, NC); // instead of uvdag(vt, wt)
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        save_vec(&up[iv2], &ut, 1);
      }
    }

  }

}

//====================================================================
template <typename REALTYPE>
void multadd_Gdn(AField<REALTYPE,SIMD2>& u, const int exu,
                 const AField<REALTYPE,SIMD2>& v, const int exv,
                 const AField<REALTYPE,SIMD2>& w, const int exw,
                 const REALTYPE a)
{
  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NDF, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NDF, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, ut2, vt[NC], wt[NDF2];

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&wt[ic2+NC*ic1], &wp[iv2], 1);
      }
    }

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv1 = iv + VLEN * (ic2 + NC * ic1);
        load_vec(&vt[ic1], &vp[iv1], 1);
      }
      for(int ic1 = 0; ic1 < NC; ++ic1){
        mult_Vdn(ut, vt, &wt[NC*ic1], NC);
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&ut2, &up[iv2], 1);
        axpy_vec(&ut2, a, &ut, 1);
        save_vec(&up[iv2], &ut2, 1);
      }
    }

  }

}

//====================================================================
template <typename REALTYPE>
void mult_Gdn(AField<REALTYPE,SIMD2>& u, const int exu,
              const AField<REALTYPE,SIMD2>& v, const int exv,
              const AField<REALTYPE,SIMD2>& w, const int exw)
{
  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NDF, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NDF, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, vt[NC], wt[NDF2];

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&wt[ic2+NC*ic1], &wp[iv2], 1);
      }
    }

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv1 = iv + VLEN * (ic2 + NC * ic1);
        load_vec(&vt[ic1], &vp[iv1], 1);
      }
      for(int ic1 = 0; ic1 < NC; ++ic1){
        mult_Vdn(ut, vt, &wt[NC*ic1], NC);
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        save_vec(&up[iv2], &ut, 1);
      }
    }

  }

}

//====================================================================
template <typename REALTYPE>
void mult_Gdd(AField<REALTYPE,SIMD2>& u, const int exu,
              const AField<REALTYPE,SIMD2>& v, const int exv,
              const AField<REALTYPE,SIMD2>& w, const int exw)
{
  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst = v.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, exu));
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(index.idx(0, NDF, 0, exv));
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(index.idx(0, NDF, 0, exw));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut, vt[NC], wt[NDF2];

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        load_vec(&wt[ic1+NC*ic2], &wp[iv2], 1);
      }
    }
    conjg_vec(wt, NDF2);

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = 0; ic1 < NC; ++ic1){
        int iv1 = iv + VLEN * (ic2 + NC * ic1);
        load_vec(&vt[ic1], &vp[iv1], 1);
      }
      for(int ic1 = 0; ic1 < NC; ++ic1){
        mult_Vdn(ut, vt, &wt[NC*ic1], NC);
        int iv2 = iv + VLEN * (ic1 + NC * ic2);
        save_vec(&up[iv2], &ut, 1);
      }
    }

  }

}

//====================================================================
// anti-hermitian
template <typename REALTYPE>
void ah_G(AField<REALTYPE,SIMD2>& u, const int ex)
{
  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst  = u.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, ex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut1, ut2;

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    for(int ic1 = 0; ic1 < NC; ++ic1){
      int iv1 = iv + VLEN * (ic1 + NC * ic1);
      load_vec(&ut1, &up[iv1], 1);
      imag_vec(&ut1, 1);
      save_vec(&up[iv1], &ut1, 1);
    }

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = ic2+1; ic1 < NC; ++ic1){
        int iv1 = iv + VLEN * (ic1 + NC * ic2);
        int iv2 = iv + VLEN * (ic2 + NC * ic1);
        load_vec(&ut1, &up[iv1], 1);
        load_vec(&ut2, &up[iv2], 1);
        conjg_vec(&ut2, 1);
        axpy_vec(&ut1, real_t(-1.0), &ut2, 1);
        scal_vec(&ut1, real_t(0.5), 1);
        set_vec(&ut2, real_t(-1.0), &ut1, 1);
        conjg_vec(&ut2, 1);
        save_vec(&up[iv1], &ut1, 1);
        save_vec(&up[iv2], &ut2, 1);
      }
    }

  }

}

//====================================================================
// anti-hermitian traceless
template <typename REALTYPE>
void at_G(AField<REALTYPE,SIMD2>& u, const int ex)
{
  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst  = u.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, ex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut1, ut2;
  real_t RNc = 1.0/real_t(NC);

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    clear_vec(&ut2, 1);
    for(int ic1 = 0; ic1 < NC; ++ic1){
      int iv1 = iv + VLEN * (ic1 + NC * ic1);
      load_vec(&ut1, &up[iv1], 1);
      add_vec(&ut2, &ut1, 1);
    }
    imag_vec(&ut2, 1);

    for(int ic1 = 0; ic1 < NC; ++ic1){
      int iv1 = iv + VLEN * (ic1 + NC * ic1);
      load_vec(&ut1, &up[iv1], 1);
      imag_vec(&ut1, 1);
      axpy_vec(&ut1, -RNc, &ut2, 1);
      save_vec(&up[iv1], &ut1, 1);
    }

    for(int ic2 = 0; ic2 < NC; ++ic2){
      for(int ic1 = ic2+1; ic1 < NC; ++ic1){

        int iv1 = iv + VLEN * (ic1 + NC * ic2);
        int iv2 = iv + VLEN * (ic2 + NC * ic1);
        load_vec(&ut1, &up[iv1], 1);
        load_vec(&ut2, &up[iv2], 1);
        conjg_vec(&ut2, 1);
        axpy_vec(&ut1, real_t(-1.0), &ut2, 1);
        scal_vec(&ut1, real_t(0.5), 1);
        set_vec(&ut2, real_t(-1.0), &ut1, 1);
        conjg_vec(&ut2, 1);
        save_vec(&up[iv1], &ut1, 1);
        save_vec(&up[iv2], &ut2, 1);
      }
    }

  }

}

//====================================================================
template <typename REALTYPE>
void add_unit(AField<REALTYPE,SIMD2>& u, const int ex, REALTYPE a)
{         // u = u + a * I (I: unit matrix)

  typedef AField<REALTYPE,SIMD2> AFIELD;
  typedef REALTYPE real_t;
  int Nst  = u.nvol();
  int Nstv = Nst/VLEN2;

  Index_lex_alt<real_t,AFIELD::IMPL> index;

  real_t* up = u.ptr(index.idx(0, NDF, 0, ex));

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  Vsimd_t ut1, ut2;
  unit_vec(&ut2, 1);

  for(int site = is; site < ns; ++site){
    int iv = VLEN * NDF2 * site;

    for(int ic1 = 0; ic1 < NC; ++ic1){
      int iv1 = iv + VLEN * (ic1 + NC * ic1);
      load_vec(&ut1, &up[iv1], 1);
      axpy_vec(&ut1, a, &ut2, 1);
      save_vec(&up[iv1], &ut1, 1);
    }

  }

}

//============================================================END=====
#endif
