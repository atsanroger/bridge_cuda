/*!
      @file    afopr_Wilson_vector-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

// This code explicitly assumes SU(3) gauge group.
#include "lib_alt_Vector/inline/define_params.h"
#include "lib_alt_Vector/inline/define_index.h"

#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)

//====================================================================
void mult_wilson_gm5_dirac(real_t *RESTRICT v2, real_t *RESTRICT v1,
                           int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_SP_R(Nst, ic, 0, ist)] = v1[IDXV_SP_R(Nst, ic, 2, ist)];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] = v1[IDXV_SP_I(Nst, ic, 2, ist)];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] = v1[IDXV_SP_R(Nst, ic, 3, ist)];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] = v1[IDXV_SP_I(Nst, ic, 3, ist)];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] = v1[IDXV_SP_R(Nst, ic, 0, ist)];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] = v1[IDXV_SP_I(Nst, ic, 0, ist)];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] = v1[IDXV_SP_R(Nst, ic, 1, ist)];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] = v1[IDXV_SP_I(Nst, ic, 1, ist)];
    }
  }

}

//====================================================================
void mult_wilson_gm5_aypx_dirac(real_t a, real_t *RESTRICT v2,
                                real_t *RESTRICT v1, int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    for(int ic = 0; ic < NC; ++ic){
      real_t v2r = v2[IDXV_SP_R(Nst, ic, 2, ist)];
      real_t v2i = v2[IDXV_SP_I(Nst, ic, 2, ist)];
      real_t v3r = v2[IDXV_SP_R(Nst, ic, 3, ist)];
      real_t v3i = v2[IDXV_SP_I(Nst, ic, 3, ist)];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] = a * v2[IDXV_SP_R(Nst, ic, 0, ist)]
                                    + v1[IDXV_SP_R(Nst, ic, 0, ist)];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] = a * v2[IDXV_SP_I(Nst, ic, 0, ist)]
                                    + v1[IDXV_SP_I(Nst, ic, 0, ist)];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] = a * v2[IDXV_SP_R(Nst, ic, 1, ist)]
                                    + v1[IDXV_SP_R(Nst, ic, 1, ist)];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] = a * v2[IDXV_SP_I(Nst, ic, 1, ist)]
                                    + v1[IDXV_SP_I(Nst, ic, 1, ist)];
      v2[IDXV_SP_R(Nst, ic, 0, ist)] = a * v2r + v1[IDXV_SP_R(Nst, ic, 2, ist)];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] = a * v2i + v1[IDXV_SP_I(Nst, ic, 2, ist)];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] = a * v3r + v1[IDXV_SP_R(Nst, ic, 3, ist)];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] = a * v3i + v1[IDXV_SP_I(Nst, ic, 3, ist)];
   }
  }

}

//====================================================================
void mult_wilson_gm5_chiral(real_t *RESTRICT v2, real_t *RESTRICT v1,
                            int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_SP_R(Nst, ic, 0, ist)] =  v1[IDXV_SP_R(Nst, ic, 0, ist)];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] =  v1[IDXV_SP_I(Nst, ic, 0, ist)];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] =  v1[IDXV_SP_R(Nst, ic, 1, ist)];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] =  v1[IDXV_SP_I(Nst, ic, 1, ist)];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] = -v1[IDXV_SP_R(Nst, ic, 2, ist)];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] = -v1[IDXV_SP_I(Nst, ic, 2, ist)];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] = -v1[IDXV_SP_R(Nst, ic, 3, ist)];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] = -v1[IDXV_SP_I(Nst, ic, 3, ist)];
    }
  }

}

//====================================================================
void mult_wilson_gm5_aypx_chiral(real_t a, real_t *RESTRICT v2,
                                 real_t *RESTRICT v1, int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    for(int id = 0; id < 2; ++id){
     for(int ic = 0; ic < NC; ++ic){
       v2[IDXV_SP_R(Nst, ic, id, ist)] =  a * v2[IDXV_SP_R(Nst, ic, id, ist)]
                                       + v1[IDXV_SP_R(Nst, ic, id, ist)];
       v2[IDXV_SP_I(Nst, ic, id, ist)] =  a * v2[IDXV_SP_I(Nst, ic, id, ist)]
                                       + v1[IDXV_SP_I(Nst, ic, id, ist)];
     }
    }
    for(int id = 2; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       v2[IDXV_SP_R(Nst, ic, id, ist)] = -a * v2[IDXV_SP_R(Nst, ic, id, ist)]
                                       - v1[IDXV_SP_R(Nst, ic, id, ist)];
       v2[IDXV_SP_I(Nst, ic, id, ist)] = -a * v2[IDXV_SP_I(Nst, ic, id, ist)]
                                       - v1[IDXV_SP_I(Nst, ic, id, ist)];
     }
    }
  }

}

//====================================================================
void mult_wilson_aypx(real_t a, real_t *RESTRICT v2, real_t *RESTRICT v1,
                      int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
   for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_SP_R(Nst, ic, id, ist)] = a * v2[IDXV_SP_R(Nst, ic, id, ist)]
                                     + v1[IDXV_SP_R(Nst, ic, id, ist)];
      v2[IDXV_SP_I(Nst, ic, id, ist)] = a * v2[IDXV_SP_I(Nst, ic, id, ist)]
                                     + v1[IDXV_SP_I(Nst, ic, id, ist)];
    }
   }
  }

}

//====================================================================
void mult_wilson_axpy(real_t *RESTRICT v2, real_t a, real_t *RESTRICT v1,
                      int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
   for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      v2[IDXV_SP_R(Nst, ic, id, ist)] += a * v1[IDXV_SP_R(Nst, ic, id, ist)];
      v2[IDXV_SP_I(Nst, ic, id, ist)] += a * v1[IDXV_SP_I(Nst, ic, id, ist)];
    }
   }
  }

}

//====================================================================
void mult_wilson_scal(real_t *RESTRICT v, real_t a,
                      int *Nsize, int Nc)
{
  int Nst  = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < is; ++ist){
   for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      v[IDXV_SP_R(Nst, ic, id, ist)] *= a;
      v[IDXV_SP_I(Nst, ic, id, ist)] *= a;
    }
   }
  }

}

//====================================================================
void mult_wilson_clear(real_t *RESTRICT v2, int *Nsize, int Nc)
{
  int Nst = Nsize[0] * Nsize[1] * Nsize[2] * Nsize[3];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  for(int ist = 0; ist < Nst; ++ist){
    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       v2[IDXV_SP_R(Nst, ic, id, ist)] = 0.0;
       v2[IDXV_SP_I(Nst, ic, id, ist)] = 0.0;
     }
    }
  }

}

//====================================================================
void mult_wilson_xp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  //  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst = Nx * Nyzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nyzt);

  real_t bc2 = bc[0];

  //  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
  for(int iyzt = is; iyzt < ns; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;

    real_t vt[NVC*ND];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic   + NVC*id] = v1[IDXV_SP_R(Nst, ic, id, ist)];
       vt[2*ic+1 + NVC*id] = v1[IDXV_SP_I(Nst, ic, id, ist)];
     }
    }

    real_t *vt1 = &buf[      NVC*2*iyzt];
    real_t *vt2 = &buf[NVC + NVC*2*iyzt];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = bc2 *( vt[icr + ID1] - vt[ici + ID4] );
      vt1[ici] = bc2 *( vt[ici + ID1] + vt[icr + ID4] );
      vt2[icr] = bc2 *( vt[icr + ID2] - vt[ici + ID3] );
      vt2[ici] = bc2 *( vt[ici + ID2] + vt[icr + ID3] );
    }

  }

}

//====================================================================
void mult_wilson_xp2(real_t *RESTRICT v2, real_t *RESTRICT u,
                     real_t *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nyzt);

  //  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
  for(int iyzt = is; iyzt < ns; ++iyzt){
    int ix  = Nx-1;
    int ist = ix + Nx * iyzt;

    real_t vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = buf[ivc +       2*NVC*iyzt];
      vt2[ivc] = buf[ivc + NVC + 2*NVC*iyzt];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] +=  wt2[1];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += -wt2[0];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] +=  wt1[1];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] += -wt1[0];

    }

  }

}

//====================================================================
 void mult_wilson_xpb(real_t *RESTRICT v2, real_t *RESTRICT u,
                      real_t *RESTRICT v1,
                      int *Nsize, int *bc, int Nc)
{
  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    int ix   = ist % Nx;
    int iyzt = ist/Nx;
    int nei  = ((ix+1) % Nx) + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == Nx-1) bc2 = bc[0];

    real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = bc2 *( vt[icr + ID1] - vt[ici + ID4] );
      vt1[ici] = bc2 *( vt[ici + ID1] + vt[icr + ID4] );
      vt2[icr] = bc2 *( vt[icr + ID2] - vt[ici + ID3] );
      vt2[ici] = bc2 *( vt[ici + ID2] + vt[icr + ID3] );
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] +=  wt2[1];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += -wt2[0];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] +=  wt1[1];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] += -wt1[0];

    }

  }

}

//====================================================================
void mult_wilson_xm1(real_t *RESTRICT buf, real_t *RESTRICT u,
                     real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nyzt);

  //  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
  for(int iyzt = is; iyzt < ns; ++iyzt){
    int ix  = Nx-1;
    int ist = ix + Nx * iyzt;

    real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];


    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic   + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, ist) ];
       vt[2*ic+1 + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, ist) ];
     }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = vt[icr + ID1] + vt[ici + ID4];
      vt1[ici] = vt[ici + ID1] - vt[icr + ID4];
      vt2[icr] = vt[icr + ID2] + vt[ici + ID3];
      vt2[ici] = vt[ici + ID2] - vt[icr + ID3];
    }

    real_t *wt1 = &buf[      NVC*2*iyzt];
    real_t *wt2 = &buf[NVC + NVC*2*iyzt];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, 0)];
      }

      wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    }

  }

}

//====================================================================
void mult_wilson_xm2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nyzt);

  real_t bc2 = bc[0];

  //  for(int iyzt = 0; iyzt < Nyzt; ++iyzt){
  for(int iyzt = is; iyzt < ns; ++iyzt){
    int ix  = 0;
    int ist = ix + Nx * iyzt;

    real_t wt1[2], wt2[2];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2 * ic;
      int ici = 2 * ic + 1;

      wt1[0] = bc2 * buf[icr       + 2*NVC*iyzt];
      wt1[1] = bc2 * buf[ici       + 2*NVC*iyzt];
      wt2[0] = bc2 * buf[icr + NVC + 2*NVC*iyzt];
      wt2[1] = bc2 * buf[ici + NVC + 2*NVC*iyzt];

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] += -wt2[1];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  wt2[0];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] += -wt1[1];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] +=  wt1[0];

    }

  }

}

//====================================================================
 void mult_wilson_xmb(real_t *RESTRICT v2, real_t *RESTRICT u,
                      real_t *RESTRICT v1,
                      int *Nsize, int *bc, int Nc)
{
  int idir = 0;

  int Nx   = Nsize[0];
  int Nyzt = Nsize[1] * Nsize[2] * Nsize[3];
  int Nst  = Nx * Nyzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    int ix   = ist % Nx;
    int iyzt = ist/Nx;
    int nei  = ix-1 + Nx * iyzt;
    if(ix == 0) nei = Nx-1 + Nx * iyzt;
    real_t bc2 = 1.0;
    if(ix == 0) bc2 = bc[0];

    real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = vt[icr + ID1] + vt[ici + ID4];
      vt1[ici] = vt[ici + ID1] - vt[icr + ID4];
      vt2[icr] = vt[icr + ID2] + vt[ici + ID3];
      vt2[ici] = vt[ici + ID2] - vt[icr + ID3];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, nei, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, nei, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  bc2 * wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  bc2 * wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  bc2 * wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  bc2 * wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] += -bc2 * wt2[1];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  bc2 * wt2[0];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] += -bc2 * wt1[1];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] +=  bc2 * wt1[0];

    }

  }

}

//====================================================================
void mult_wilson_yp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nx*Nzt);

  real_t bc2 = bc[1];

  //  for(int izt = 0; izt < Nzt; ++izt){
  //   for(int ix = 0; ix < Nx; ++ix){
  for(int ixzt = is; ixzt < ns; ++ixzt){
   {
     int ix  = ixzt % Nx;
     int izt = ixzt/Nx;

     int iy   = 0;
     int ist  = ix + Nx * (iy + Ny * izt);
     int ixzt = ix + Nx * izt;

     real_t vt[NVC*ND];

     for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        vt[2*ic   + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, ist) ];
        vt[2*ic+1 + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, ist) ];
      }
     }

     real_t *vt1 = &buf[      NVC*2*ixzt];
     real_t *vt2 = &buf[NVC + NVC*2*ixzt];

     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = bc2 *( vt[icr + ID1] + vt[icr + ID4] );
       vt1[ici] = bc2 *( vt[ici + ID1] + vt[ici + ID4] );
       vt2[icr] = bc2 *( vt[icr + ID2] - vt[icr + ID3] );
       vt2[ici] = bc2 *( vt[ici + ID2] - vt[ici + ID3] );
     }

   }
  }

}

//====================================================================
void mult_wilson_yp2(real_t *RESTRICT v2, real_t *RESTRICT u,
                     real_t *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nx*Nzt);

  //  for(int izt = 0; izt < Nzt; ++izt){
  //   for(int ix = 0; ix < Nx; ++ix){
  for(int ixzt = is; ixzt < ns; ++ixzt){
   {
     int ix  = ixzt % Nx;
     int izt = ixzt/Nx;

     int iy   = Ny-1;
     int ist  = ix + Nx * (iy + Ny * izt);
     int ixzt = ix + Nx * izt;

     real_t vt1[NVC], vt2[NVC], ut[NDF];
     real_t wt1[2], wt2[2];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = buf[ivc +       2*NVC*ixzt];
      vt2[ivc] = buf[ivc + NVC + 2*NVC*ixzt];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] += -wt2[0];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += -wt2[1];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] +=  wt1[1];

    }

   }
  }

}

//====================================================================
 void mult_wilson_ypb(real_t *RESTRICT v2, real_t *RESTRICT u,
                      real_t *RESTRICT v1,
                      int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  for(int ist = is; ist < ns; ++ist){
    int ix  = ist % Nx;
    int iy  = (ist/Nx) % Ny;
    int izt = ist/(Nx*Ny);
    int iy2 = (iy+1) % Ny;
    int nei = ix + Nx * (iy2 + Ny * izt);
    real_t bc2 = 1.0;
    if(iy == Ny-1) bc2 = bc[idir];

    real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];


    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
       vt1[icr] = bc2 *( vt[icr + ID1] + vt[icr + ID4] );
       vt1[ici] = bc2 *( vt[ici + ID1] + vt[ici + ID4] );
       vt2[icr] = bc2 *( vt[icr + ID2] - vt[icr + ID3] );
       vt2[ici] = bc2 *( vt[ici + ID2] - vt[ici + ID3] );
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] += -wt2[0];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += -wt2[1];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] +=  wt1[1];
    }

  }

}

//====================================================================
void mult_wilson_ym1(real_t *RESTRICT buf, real_t *RESTRICT u,
                     real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nzt = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nx*Nzt);

  //  for(int izt = 0; izt < Nzt; ++izt){
  //   for(int ix = 0; ix < Nx; ++ix){
  for(int ixzt = is; ixzt < ns; ++ixzt){
   {
     int ix  = ixzt % Nx;
     int izt = ixzt/Nx;

     int iy   = Ny-1;
     int ist  = ix + Nx * (iy + Ny*izt);
     int ixzt = ix + Nx * izt;

     real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];

     for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        vt[2*ic   + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, ist) ];
        vt[2*ic+1 + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, ist) ];
      }
     }

     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = vt[icr + ID1] - vt[icr + ID4];
       vt1[ici] = vt[ici + ID1] - vt[ici + ID4];
       vt2[icr] = vt[icr + ID2] + vt[icr + ID3];
       vt2[ici] = vt[ici + ID2] + vt[ici + ID3];
     }

     real_t *wt1 = &buf[      NVC*2*ixzt];
     real_t *wt2 = &buf[NVC + NVC*2*ixzt];

     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;

       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, 0)];
         ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, 0)];
       }

       wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt2[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
       wt2[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
     }

   }
  }

}

//====================================================================
void mult_wilson_ym2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nst = Nx * Ny * Nzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nx*Nzt);

  real_t bc2 = bc[idir];

  //  for(int izt = 0; izt < Nzt; ++izt){
  //   for(int ix = 0; ix < Nx; ++ix){
  for(int ixzt = is; ixzt < ns; ++ixzt){
   {
     int ix  = ixzt % Nx;
     int izt = ixzt/Nx;

     int iy = 0;
     int ist  = ix + Nx*(iy + Ny*izt);
     int ixzt = ix + Nx * izt;

     real_t wt1[2], wt2[2];

     for(int ic = 0; ic < NC; ++ic){
       int icr = 2 * ic;
       int ici = 2 * ic + 1;

       wt1[0] = bc2 * buf[icr       + 2*NVC*ixzt];
       wt1[1] = bc2 * buf[ici       + 2*NVC*ixzt];
       wt2[0] = bc2 * buf[icr + NVC + 2*NVC*ixzt];
       wt2[1] = bc2 * buf[ici + NVC + 2*NVC*ixzt];

       v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
       v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
       v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
       v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
       v2[IDXV_SP_R(Nst, ic, 2, ist)] +=  wt2[0];
       v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  wt2[1];
       v2[IDXV_SP_R(Nst, ic, 3, ist)] += -wt1[0];
       v2[IDXV_SP_I(Nst, ic, 3, ist)] += -wt1[1];
    }

   }
  }

}

//====================================================================
 void mult_wilson_ymb(real_t *RESTRICT v2, real_t *RESTRICT u,
                      real_t *RESTRICT v1,
                      int *Nsize, int *bc, int Nc)
{
  int idir = 1;

  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nzt  = Nsize[2] * Nsize[3];
  int Nst  = Nx * Ny * Nzt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    int ix  = ist % Nx;
    int iy  = (ist/Nx) % Ny;
    int izt = ist/(Nx*Ny);
    int iy2 = (iy-1+Ny) % Ny;
    int nei = ix + Nx * (iy2 + Ny * izt);

    real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

     real_t bc2 = 1.0;
     if(iy == 0) bc2 = bc[idir];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = vt[icr + ID1] - vt[icr + ID4];
      vt1[ici] = vt[ici + ID1] - vt[ici + ID4];
      vt2[icr] = vt[icr + ID2] + vt[icr + ID3];
      vt2[ici] = vt[ici + ID2] + vt[ici + ID3];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, nei, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, nei, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  bc2 * wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  bc2 * wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  bc2 * wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  bc2 * wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] +=  bc2 * wt2[0];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  bc2 * wt2[1];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] += -bc2 * wt1[0];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] += -bc2 * wt1[1];
    }

  }

}

//====================================================================
void mult_wilson_zp1(real_t *RESTRICT buf, real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxy*Nt);

  real_t bc2 = bc[idir];

  //  for(int it = 0; it < Nt; ++it){
  //   for(int ixy = 0; ixy < Nxy; ++ixy){
 for(int ixyt = is; ixyt < ns; ++ixyt){
   {
     int ixy = ixyt % Nxy;
     int it  = ixyt/Nxy;

     int iz   = 0;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     real_t vt[NVC*ND];

     for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        vt[2*ic   + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, ist) ];
        vt[2*ic+1 + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, ist) ];
      }
     }

     real_t *vt1 = &buf[      NVC*2*ixyt];
     real_t *vt2 = &buf[NVC + NVC*2*ixyt];

     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = bc2 *( vt[icr + ID1] - vt[ici + ID3] );
       vt1[ici] = bc2 *( vt[ici + ID1] + vt[icr + ID3] );
       vt2[icr] = bc2 *( vt[icr + ID2] + vt[ici + ID4] );
       vt2[ici] = bc2 *( vt[ici + ID2] - vt[icr + ID4] );
     }

   }
  }

}

//====================================================================
void mult_wilson_zp2(real_t *RESTRICT v2, real_t *RESTRICT u,
                     real_t *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxy*Nt);

  //  for(int it = 0; it < Nt; ++it){
  //   for(int ixy = 0; ixy < Nxy; ++ixy){
 for(int ixyt = is; ixyt < ns; ++ixyt){
   {
     int ixy = ixyt % Nxy;
     int it  = ixyt/Nxy;

     int iz = Nz-1;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     real_t vt1[NVC], vt2[NVC], ut[NDF];
     real_t wt1[2], wt2[2];


    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = buf[ivc +       2*NVC*ixyt];
      vt2[ivc] = buf[ivc + NVC + 2*NVC*ixyt];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] +=  wt1[1];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += -wt1[0];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] += -wt2[1];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] +=  wt2[0];

    }

   }
  }

}

//====================================================================
 void mult_wilson_zpb(real_t *RESTRICT v2, real_t *RESTRICT u,
                      real_t *RESTRICT v1,
                      int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  for(int ist = is; ist < ns; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz+1) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);
    real_t bc2 = 1.0;
    if(iz == Nz-1) bc2 = bc[idir];

    real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
       vt1[icr] = bc2 *( vt[icr + ID1] - vt[ici + ID3] );
       vt1[ici] = bc2 *( vt[ici + ID1] + vt[icr + ID3] );
       vt2[icr] = bc2 *( vt[icr + ID2] + vt[ici + ID4] );
       vt2[ici] = bc2 *( vt[ici + ID2] - vt[icr + ID4] );
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] +=  wt1[1];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += -wt1[0];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] += -wt2[1];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] +=  wt2[0];
    }

  }

}

//====================================================================
void mult_wilson_zm1(real_t *RESTRICT buf, real_t *RESTRICT u,
                     real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxy*Nt);

  //  for(int it = 0; it < Nt; ++it){
  //   for(int ixy = 0; ixy < Nxy; ++ixy){
 for(int ixyt = is; ixyt < ns; ++ixyt){
   {
     int ixy = ixyt % Nxy;
     int it  = ixyt/Nxy;

     int iz = Nz-1;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];

     for(int id = 0; id < ND; ++id){
      for(int ic = 0; ic < NC; ++ic){
        vt[2*ic   + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, ist) ];
        vt[2*ic+1 + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, ist) ];
      }
     }

     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt1[icr] = vt[icr + ID1] + vt[ici + ID3];
       vt1[ici] = vt[ici + ID1] - vt[icr + ID3];
       vt2[icr] = vt[icr + ID2] - vt[ici + ID4];
       vt2[ici] = vt[ici + ID2] + vt[icr + ID4];
     }

     real_t *wt1 = &buf[      NVC*2*ixyt];
     real_t *wt2 = &buf[NVC + NVC*2*ixyt];

     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;

       for(int ic2 = 0; ic2 < NC; ++ic2){
         ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, 0)];
         ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, 0)];
       }

       wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
       wt2[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
       wt2[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                            vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    }

   }
  }

}

//====================================================================
void mult_wilson_zm2(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  real_t bc2 = bc[idir];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxy*Nt);

  //  for(int it = 0; it < Nt; ++it){
  //   for(int ixy = 0; ixy < Nxy; ++ixy){
 for(int ixyt = is; ixyt < ns; ++ixyt){
   {
     int ixy = ixyt % Nxy;
     int it  = ixyt/Nxy;

     int iz = 0;
     int ist  = ixy + Nxy * (iz + Nz * it);
     int ixyt = ixy + Nxy * it;

     real_t wt1[2], wt2[2];

     for(int ic = 0; ic < NC; ++ic){
       int icr = 2 * ic;
       int ici = 2 * ic + 1;

       wt1[0] = bc2 * buf[icr       + 2*NVC*ixyt];
       wt1[1] = bc2 * buf[ici       + 2*NVC*ixyt];
       wt2[0] = bc2 * buf[icr + NVC + 2*NVC*ixyt];
       wt2[1] = bc2 * buf[ici + NVC + 2*NVC*ixyt];

       v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
       v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
       v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
       v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
       v2[IDXV_SP_R(Nst, ic, 2, ist)] += -wt1[1];
       v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  wt1[0];
       v2[IDXV_SP_R(Nst, ic, 3, ist)] +=  wt2[1];
       v2[IDXV_SP_I(Nst, ic, 3, ist)] += -wt2[0];
     }

   }
  }

}

//==================================================================== 
void mult_wilson_zmb(real_t *RESTRICT v2, real_t *RESTRICT u,
                     real_t *RESTRICT v1,
                     int *Nsize, int *bc, int Nc)
{
  int idir = 2;

  int Nxy  = Nsize[0] * Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxy * Nz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    int ixy = ist % Nxy;
    int iz  = (ist/Nxy) % Nz;
    int it  = ist/(Nxy*Nz);
    int iz2 = (iz-1+Nz) % Nz;
    int nei = ixy + Nxy * (iz2 + Nz * it);
    real_t bc2 = 1.0;
    if(iz == 0) bc2 = bc[idir];

    real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];


    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = vt[icr + ID1] + vt[ici + ID3];
      vt1[ici] = vt[ici + ID1] - vt[icr + ID3];
      vt2[icr] = vt[icr + ID2] - vt[ici + ID4];
      vt2[ici] = vt[ici + ID2] + vt[icr + ID4];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, nei, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, nei, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  bc2 * wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  bc2 * wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  bc2 * wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  bc2 * wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] += -bc2 * wt1[1];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] +=  bc2 * wt1[0];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] +=  bc2 * wt2[1];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] += -bc2 * wt2[0];
    }

  }

}

//====================================================================
void mult_wilson_tp1_dirac(real_t *RESTRICT buf, real_t *RESTRICT v1,
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  real_t bc2 = bc[idir];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxyz);

  //  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it = 0;
    int ist  = ixyz + Nxyz * it;

    real_t vt[NVC*ND];

    //for(int id = 0; id < ND; ++id){
    for(int id = 2; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic   + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, ist) ];
       vt[2*ic+1 + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, ist) ];
     }
    }

    real_t *vt1 = &buf[      NVC*2*ixyz];
    real_t *vt2 = &buf[NVC + NVC*2*ixyz];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = bc2 * 2.0 * vt[icr + ID3];
      vt1[ici] = bc2 * 2.0 * vt[ici + ID3];
      vt2[icr] = bc2 * 2.0 * vt[icr + ID4];
      vt2[ici] = bc2 * 2.0 * vt[ici + ID4];
    }

  }

}

//====================================================================
void mult_wilson_tp2_dirac(real_t *RESTRICT v2, real_t *RESTRICT u,
                           real_t *RESTRICT buf, 
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxyz);

  //  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

    real_t vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = buf[ivc +       2*NVC*ixyz];
      vt2[ivc] = buf[ivc + NVC + 2*NVC*ixyz];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 2, ist)] += wt1[0];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += wt1[1];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] += wt2[0];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] += wt2[1];

    }

  }

}

//====================================================================
 void mult_wilson_tpb_dirac(real_t *RESTRICT v2, real_t *RESTRICT u,
                            real_t *RESTRICT v1,
                            int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    int ixyz = ist % Nxyz;
    int it  = ist/Nxyz;
    int it2 = (it+1) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == Nt-1) bc2 = bc[idir];

    real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];


    //for(int id = 0; id < ND; ++id){
    for(int id = 2; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
       vt1[icr] = bc2 * 2.0 * vt[icr + ID3];
       vt1[ici] = bc2 * 2.0 * vt[ici + ID3];
       vt2[icr] = bc2 * 2.0 * vt[icr + ID4];
       vt2[ici] = bc2 * 2.0 * vt[ici + ID4];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 2, ist)] += wt1[0];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += wt1[1];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] += wt2[0];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] += wt2[1];
    }

  }

}

//====================================================================
void mult_wilson_tm1_dirac(real_t *RESTRICT buf, real_t *RESTRICT u,
                           real_t *RESTRICT v1,
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxyz);

  //  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

    real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];

    //for(int id = 0; id < ND; ++id){
    for(int id = 0; id < 2; ++id){
     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic   + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, ist) ];
       vt[2*ic+1 + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, ist) ];
     }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = 2.0 * vt[icr + ID1];
      vt1[ici] = 2.0 * vt[ici + ID1];
      vt2[icr] = 2.0 * vt[icr + ID2];
      vt2[ici] = 2.0 * vt[ici + ID2];
    }

    real_t *wt1 = &buf[      NVC*2*ixyz];
    real_t *wt2 = &buf[NVC + NVC*2*ixyz];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, 0)];
      }

      wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    }

  }

}

//====================================================================
void mult_wilson_tm2_dirac(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxyz);

  real_t bc2 = bc[idir];

  //  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it  = 0;
    int ist = ixyz + Nxyz * it;

    real_t wt1[2], wt2[2];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2 * ic;
      int ici = 2 * ic + 1;

      wt1[0] = bc2 * buf[icr       + 2 * NVC * ixyz];
      wt1[1] = bc2 * buf[ici       + 2 * NVC * ixyz];
      wt2[0] = bc2 * buf[icr + NVC + 2 * NVC * ixyz];
      wt2[1] = bc2 * buf[ici + NVC + 2 * NVC * ixyz];

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
    }

  }

}

//==================================================================== 
void mult_wilson_tmb_dirac(real_t *RESTRICT v2, real_t *RESTRICT u,
                           real_t *RESTRICT v1,
                           int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    int ixyz = ist % Nxyz;
    int it   = ist/Nxyz;
    int it2 = (it-1+Nt) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == 0) bc2 = bc[idir];

    real_t vt[2*NC*ND];
    real_t vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    //for(int id = 0; id < ND; ++id){
    for(int id = 0; id < 2; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt1[icr] = 2.0 * vt[icr + ID1];
      vt1[ici] = 2.0 * vt[ici + ID1];
      vt2[icr] = 2.0 * vt[icr + ID2];
      vt2[ici] = 2.0 * vt[ici + ID2];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, nei, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, nei, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  bc2 * wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  bc2 * wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  bc2 * wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  bc2 * wt2[1];
    }

  }

}

//====================================================================
void mult_wilson_tp1_chiral(real_t *RESTRICT buf, real_t *RESTRICT v1,
                            int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  real_t bc2 = bc[idir];

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxyz);

  //  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it = 0;
    int ist  = ixyz + Nxyz * it;

    real_t vt[NVC*ND];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic   + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, ist) ];
       vt[2*ic+1 + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, ist) ];
     }
    }

    real_t *vt1 = &buf[      NVC*2*ixyz];
    real_t *vt2 = &buf[NVC + NVC*2*ixyz];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = bc2 * (vt[ivc + ID1] + vt[ivc + ID3]);
      vt2[ivc] = bc2 * (vt[ivc + ID2] + vt[ivc + ID4]);
    }

  }

}

//====================================================================
void mult_wilson_tp2_chiral(real_t *RESTRICT v2, real_t *RESTRICT u,
                            real_t *RESTRICT buf, 
                            int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxyz);

  //  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

    real_t vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = buf[ivc +       2*NVC*ixyz];
      vt2[ivc] = buf[ivc + NVC + 2*NVC*ixyz];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] += wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] += wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] += wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] += wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] += wt1[0];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += wt1[1];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] += wt2[0];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] += wt2[1];

    }

  }

}

//====================================================================
 void mult_wilson_tpb_chiral(real_t *RESTRICT v2, real_t *RESTRICT u,
                             real_t *RESTRICT v1,
                             int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    int ixyz = ist % Nxyz;
    int it  = ist/Nxyz;
    int it2 = (it+1) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == Nt-1) bc2 = bc[idir];

    real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];


    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = bc2 * (vt[ivc + ID1] + vt[ivc + ID3]);
      vt2[ivc] = bc2 * (vt[ivc + ID2] + vt[ivc + ID4]);
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] = u[IDXV_G_R(Nst, ic2, ic, ist, 0)];
        ut[2*ic2+1] = u[IDXV_G_I(Nst, ic2, ic, ist, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] += wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] += wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] += wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] += wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] += wt1[0];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] += wt1[1];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] += wt2[0];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] += wt2[1];
    }

  }

}

//====================================================================
void mult_wilson_tm1_chiral(real_t *RESTRICT buf, real_t *RESTRICT u,
                            real_t *RESTRICT v1,
                            int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxyz);

  //  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it   = Nt-1;
    int ist  = ixyz + Nxyz * it;

    real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       vt[2*ic   + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, ist) ];
       vt[2*ic+1 + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, ist) ];
     }
    }

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = vt[ivc + ID1] - vt[ivc + ID3];
      vt2[ivc] = vt[ivc + ID2] - vt[ivc + ID4];
    }

    real_t *wt1 = &buf[      NVC*2*ixyz];
    real_t *wt2 = &buf[NVC + NVC*2*ixyz];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, ist, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, ist, 0)];
      }

      wt1[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[icr] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[ici] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                           vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    }

  }

}

//====================================================================
void mult_wilson_tm2_chiral(real_t *RESTRICT v2, real_t *RESTRICT buf, 
                            int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nxyz);

  real_t bc2 = bc[idir];

  //  for(int ixyz = 0; ixyz < Nxyz; ++ixyz){
  for(int ixyz = is; ixyz < ns; ++ixyz){
    int it  = 0;
    int ist = ixyz + Nxyz * it;

    real_t wt1[2], wt2[2];

    for(int ic = 0; ic < NC; ++ic){
      int icr = 2 * ic;
      int ici = 2 * ic + 1;

      wt1[0] = bc2 * buf[icr       + 2 * NVC * ixyz];
      wt1[1] = bc2 * buf[ici       + 2 * NVC * ixyz];
      wt2[0] = bc2 * buf[icr + NVC + 2 * NVC * ixyz];
      wt2[1] = bc2 * buf[ici + NVC + 2 * NVC * ixyz];

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] -=  wt1[0];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] -=  wt1[1];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] -=  wt2[0];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] -=  wt2[1];
    }

  }

}

//==================================================================== 
void mult_wilson_tmb_chiral(real_t *RESTRICT v2, real_t *RESTRICT u,
                            real_t *RESTRICT v1,
                            int *Nsize, int *bc, int Nc)
{
  int idir = 3;

  int Nxyz = Nsize[0] * Nsize[1] * Nsize[2];
  int Nt   = Nsize[3];
  int Nst = Nxyz * Nt;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, Nst);

  //  for(int ist = 0; ist < Nst; ++ist){
  for(int ist = is; ist < ns; ++ist){
    int ixyz = ist % Nxyz;
    int it   = ist/Nxyz;
    int it2 = (it-1+Nt) % Nt;
    int nei = ixyz + Nxyz * it2;
    real_t bc2 = 1.0;
    if(it == 0) bc2 = bc[idir];

    real_t vt[2*NC*ND];
    real_t vt1[NVC], vt2[NVC], ut[NDF];
    real_t wt1[2], wt2[2];

    for(int id = 0; id < ND; ++id){
     for(int ic = 0; ic < NC; ++ic){
       int icr = 2*ic;
       int ici = 2*ic + 1;
       vt[icr + NVC*id] = v1[ IDXV_SP_R(Nst, ic, id, nei) ];
       vt[ici + NVC*id] = v1[ IDXV_SP_I(Nst, ic, id, nei) ];
      }
    }

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt1[ivc] = vt[ivc + ID1] - vt[ivc + ID3];
      vt2[ivc] = vt[ivc + ID2] - vt[ivc + ID4];
    }

    for(int ic = 0; ic < NC; ++ic){

      for(int ic2 = 0; ic2 < NC; ++ic2){
        ut[2*ic2  ] =   u[IDXV_G_R(Nst, ic, ic2, nei, 0)];
        ut[2*ic2+1] = - u[IDXV_G_I(Nst, ic, ic2, nei, 0)];
      }

      wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
      wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
      wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                         vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

      v2[IDXV_SP_R(Nst, ic, 0, ist)] +=  bc2 * wt1[0];
      v2[IDXV_SP_I(Nst, ic, 0, ist)] +=  bc2 * wt1[1];
      v2[IDXV_SP_R(Nst, ic, 1, ist)] +=  bc2 * wt2[0];
      v2[IDXV_SP_I(Nst, ic, 1, ist)] +=  bc2 * wt2[1];
      v2[IDXV_SP_R(Nst, ic, 2, ist)] -=  bc2 * wt1[0];
      v2[IDXV_SP_I(Nst, ic, 2, ist)] -=  bc2 * wt1[1];
      v2[IDXV_SP_R(Nst, ic, 3, ist)] -=  bc2 * wt2[0];
      v2[IDXV_SP_I(Nst, ic, 3, ist)] -=  bc2 * wt2[1];
    }

  }

}

//====================================================================

