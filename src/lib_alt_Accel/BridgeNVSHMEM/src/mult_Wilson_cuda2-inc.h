/*!
      @file    mult_Wilson_cuda2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2538 $
*/


#define MULT_UV_R(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v0-u1*v1 + u2*v2-u3*v3 + u4*v4-u5*v5)
#define MULT_UV_I(u0,u1,u2,u3,u4,u5,v0,v1,v2,v3,v4,v5)   (u0*v1+u1*v0 + u2*v3+u3*v2 + u4*v5+u5*v4)


//====================================================================
__global__
void mult_wilson_xp1_kernel(real_t* buf, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_x)
{
  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = 0;
  int ist = ix + Nx * iyzt;

  real_t bc2 = bc_x;

  real_t vt[NVC*ND];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
      vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
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

//====================================================================
void mult_wilson_xp1(real_t* buf, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nyzt = Ny * Nz * Nt;

  int bc_x = bc[0];

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nyzt/nth;
  mult_wilson_xp1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
				      Nx, Ny, Nz, Nt, bc_x);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_xp2_kernel(real_t* v2, real_t* u, real_t* buf, 
                            int Nx, int Ny, int Nz, int Nt)
{
  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

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
      ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
      ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] +=  wt2[1];
    v2[IDX2_SP_I(ic, 2, ist)] += -wt2[0];
    v2[IDX2_SP_R(ic, 3, ist)] +=  wt1[1];
    v2[IDX2_SP_I(ic, 3, ist)] += -wt1[0];
  }

}

//====================================================================
void mult_wilson_xp2(real_t* v2,  real_t* u, real_t* buf,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nyzt = Ny * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* buf_dev = (real_t*)dev_ptr(buf);

  int nth = VECTOR_LENGTH;
  int nbl = Nyzt/nth;
  mult_wilson_xp2_kernel<<<nbl,nth>>>(v2_dev, u_dev, buf_dev,
				      Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_xpb_kernel(real_t* v2, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_x)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix   = ist % Nx;
  int iyzt = ist / Nx;
  int nei  = ((ix+1) % Nx) + Nx * iyzt;

  real_t bc2 = 1.0;
  if(ix == Nx-1) bc2 = bc_x;

  real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
      vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
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
      ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
      ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] +=  wt2[1];
    v2[IDX2_SP_I(ic, 2, ist)] += -wt2[0];
    v2[IDX2_SP_R(ic, 3, ist)] +=  wt1[1];
    v2[IDX2_SP_I(ic, 3, ist)] += -wt1[0];
  }

}

//====================================================================
void mult_wilson_xpb(real_t* v2, real_t* u, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int bc_x = bc[0];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  mult_wilson_xpb_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
				      Nx, Ny, Nz, Nt, bc_x);

  CHECK(cudaDeviceSynchronize());



}

//====================================================================
__global__
void mult_wilson_xm1_kernel(real_t* buf, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt)
{
  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = Nx-1;
  int ist = ix + Nx * iyzt;

  real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
      vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
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
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
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

//====================================================================
void mult_wilson_xm1(real_t* buf, real_t* u, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nyzt = Ny * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nyzt/nth;
  mult_wilson_xm1_kernel<<<nbl,nth>>>(buf_dev, u_dev, v1_dev,
				      Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_xm2_kernel(real_t* v2, real_t* buf, 
                            int Nx, int Ny, int Nz, int Nt, int bc_x)
{
  real_t bc2 = bc_x;

  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

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

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] += -wt2[1];
    v2[IDX2_SP_I(ic, 2, ist)] +=  wt2[0];
    v2[IDX2_SP_R(ic, 3, ist)] += -wt1[1];
    v2[IDX2_SP_I(ic, 3, ist)] +=  wt1[0];
  }

}

//====================================================================
void mult_wilson_xm2(real_t* v2, real_t* buf,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nyzt = Ny * Nz * Nt;

  int bc_x = bc[0];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);

  int nth = VECTOR_LENGTH;
  int nbl = Nyzt/nth;
  mult_wilson_xm2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
				      Nx, Ny, Nz, Nt, bc_x);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_xmb_kernel(real_t* v2, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_x)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix   = ist % Nx;
  int iyzt = ist / Nx;

  int nei  = ix-1 + Nx * iyzt;
  if(ix == 0) nei = Nx-1 + Nx * iyzt;

  real_t bc2 = 1.0;
  if(ix == 0) bc2 = bc_x;

  real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
      vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
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
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, nei)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, nei)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  bc2 * wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  bc2 * wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  bc2 * wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  bc2 * wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] += -bc2 * wt2[1];
    v2[IDX2_SP_I(ic, 2, ist)] +=  bc2 * wt2[0];
    v2[IDX2_SP_R(ic, 3, ist)] += -bc2 * wt1[1];
    v2[IDX2_SP_I(ic, 3, ist)] +=  bc2 * wt1[0];
  }

}

//====================================================================
void mult_wilson_xmb(real_t* v2, real_t* u, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int bc_x = bc[0];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  mult_wilson_xmb_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
				      Nx, Ny, Nz, Nt, bc_x);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_yp1_kernel(real_t* buf, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_y)
{
  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ixzt % Nx;
  int izt = ixzt / Nx;
  int iy   = 0;
  int ist  = ix + Nx * (iy + Ny * izt);

  real_t bc2 = bc_y;

  real_t vt[NVC*ND];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
      vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
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

//====================================================================
void mult_wilson_yp1(real_t* buf, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxzt = Nx * Nz * Nt;

  int bc_y = bc[1];

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nxzt/nth;
  mult_wilson_yp1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
				      Nx, Ny, Nz, Nt, bc_y);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_yp2_kernel(real_t* v2, real_t* u, real_t* buf, 
                            int Nx, int Ny, int Nz, int Nt)
{
  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ixzt % Nx;
  int izt = ixzt / Nx;
  int iy  = Ny-1;
  int ist = ix + Nx * (iy + Ny * izt);

  real_t vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  for(int ivc = 0; ivc < NVC; ++ivc){
    vt1[ivc] = buf[ivc +       2*NVC*ixzt];
    vt2[ivc] = buf[ivc + NVC + 2*NVC*ixzt];
  }

  for(int ic = 0; ic < NC; ++ic){

    for(int ic2 = 0; ic2 < NC; ++ic2){
      ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
      ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] += -wt2[0];
    v2[IDX2_SP_I(ic, 2, ist)] += -wt2[1];
    v2[IDX2_SP_R(ic, 3, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 3, ist)] +=  wt1[1];
  }

}

//====================================================================
void mult_wilson_yp2(real_t* v2,  real_t* u, real_t* buf,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxzt = Nx * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* buf_dev = (real_t*)dev_ptr(buf);

  int nth = VECTOR_LENGTH;
  int nbl = Nxzt/nth;
  mult_wilson_yp2_kernel<<<nbl,nth>>>(v2_dev, u_dev, buf_dev,
				      Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_ypb_kernel(real_t* v2, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_y)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ist % Nx;
  int iy  = (ist/Nx) % Ny;
  int izt = ist/(Nx*Ny);

  int iy2 = (iy+1) % Ny;
  int nei = ix + Nx * (iy2 + Ny * izt);

  real_t bc2 = 1.0;
  if(iy == Ny-1) bc2 = bc_y;

  real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  for(int id = 0; id < ND; ++id){
   for(int ic = 0; ic < NC; ++ic){
     int icr = 2*ic;
     int ici = 2*ic + 1;
     vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
     vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
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
      ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
      ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] += -wt2[0];
    v2[IDX2_SP_I(ic, 2, ist)] += -wt2[1];
    v2[IDX2_SP_R(ic, 3, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 3, ist)] +=  wt1[1];
  }

}

//====================================================================
void mult_wilson_ypb(real_t* v2, real_t* u, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int bc_y = bc[1];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  mult_wilson_ypb_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                      Nx, Ny, Nz, Nt, bc_y);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_ym1_kernel(real_t* buf, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt)
{
  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ixzt % Nx;
  int izt = ixzt / Nx;
  int iy   = Ny-1;
  int ist  = ix + Nx * (iy + Ny*izt);

  real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
      vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
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
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
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

//====================================================================
void mult_wilson_ym1(real_t* buf, real_t* u, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxzt = Nx * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nxzt/nth;
  mult_wilson_ym1_kernel<<<nbl,nth>>>(buf_dev, u_dev, v1_dev,
				      Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_ym2_kernel(real_t* v2, real_t* buf, 
                            int Nx, int Ny, int Nz, int Nt, int bc_y)
{
  real_t bc2 = bc_y;

  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ixzt % Nx;
  int izt = ixzt / Nx;
  int iy = 0;
  int ist  = ix + Nx*(iy + Ny*izt);

  real_t wt1[2], wt2[2];

  for(int ic = 0; ic < NC; ++ic){
    int icr = 2 * ic;
    int ici = 2 * ic + 1;

    wt1[0] = bc2 * buf[icr       + 2*NVC*ixzt];
    wt1[1] = bc2 * buf[ici       + 2*NVC*ixzt];
    wt2[0] = bc2 * buf[icr + NVC + 2*NVC*ixzt];
    wt2[1] = bc2 * buf[ici + NVC + 2*NVC*ixzt];

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 2, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 3, ist)] += -wt1[0];
    v2[IDX2_SP_I(ic, 3, ist)] += -wt1[1];
  }

}

//====================================================================
void mult_wilson_ym2(real_t* v2, real_t* buf,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxzt = Nx * Nz * Nt;

  int bc_y = bc[1];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);

  int nth = VECTOR_LENGTH;
  int nbl = Nxzt/nth;
  mult_wilson_ym2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
				      Nx, Ny, Nz, Nt, bc_y);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_ymb_kernel(real_t* v2, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_y)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ist % Nx;
  int iy  = (ist/Nx) % Ny;
  int izt = ist/(Nx*Ny);
  int iy2 = (iy-1+Ny) % Ny;
  int nei = ix + Nx * (iy2 + Ny * izt);

  real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  real_t bc2 = 1.0;
  if(iy == 0) bc2 = bc_y;

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
      vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
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
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, nei)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, nei)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  bc2 * wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  bc2 * wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  bc2 * wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  bc2 * wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] +=  bc2 * wt2[0];
    v2[IDX2_SP_I(ic, 2, ist)] +=  bc2 * wt2[1];
    v2[IDX2_SP_R(ic, 3, ist)] += -bc2 * wt1[0];
    v2[IDX2_SP_I(ic, 3, ist)] += -bc2 * wt1[1];
  }

}

//====================================================================
void mult_wilson_ymb(real_t* v2, real_t* u, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int bc_y = bc[1];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  mult_wilson_ymb_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                      Nx, Ny, Nz, Nt, bc_y);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_zp1_kernel(real_t* buf, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_z)
{
  int Nxy  = Nx * Ny;

  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int ixy = ixyt % Nxy;
  int it  = ixyt / Nxy;
  int iz   = 0;
  int ist  = ixy + Nxy * (iz + Nz * it);

  real_t bc2 = bc_z;

  real_t vt[NVC*ND];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
      vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
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

//====================================================================
void mult_wilson_zp1(real_t* buf, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxyt = Nx * Ny * Nt;

  int bc_z = bc[2];

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nxyt/nth;
  mult_wilson_zp1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
				      Nx, Ny, Nz, Nt, bc_z);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_zp2_kernel(real_t* v2, real_t* u, real_t* buf, 
                            int Nx, int Ny, int Nz, int Nt)
{
  int Nxy = Nx * Ny;

  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int ixy = ixyt % Nxy;
  int it  = ixyt / Nxy;
  int iz  = Nz-1;
  int ist = ixy + Nxy * (iz + Nz * it);

  real_t vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  for(int ivc = 0; ivc < NVC; ++ivc){
    vt1[ivc] = buf[ivc +       2*NVC*ixyt];
    vt2[ivc] = buf[ivc + NVC + 2*NVC*ixyt];
  }

  for(int ic = 0; ic < NC; ++ic){

    for(int ic2 = 0; ic2 < NC; ++ic2){
      ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
      ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] +=  wt1[1];
    v2[IDX2_SP_I(ic, 2, ist)] += -wt1[0];
    v2[IDX2_SP_R(ic, 3, ist)] += -wt2[1];
    v2[IDX2_SP_I(ic, 3, ist)] +=  wt2[0];
  }

}

//====================================================================
void mult_wilson_zp2(real_t* v2,  real_t* u, real_t* buf,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxyt = Nx * Ny * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* buf_dev = (real_t*)dev_ptr(buf);

  int nth = VECTOR_LENGTH;
  int nbl = Nxyt/nth;
  mult_wilson_zp2_kernel<<<nbl,nth>>>(v2_dev, u_dev, buf_dev,
				      Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_zpb_kernel(real_t* v2, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_z)
{
  int Nxy = Nx * Ny;

  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ixy = ist % Nxy;
  int iz  = (ist/Nxy) % Nz;
  int it  = ist/(Nxy*Nz);

  int iz2 = (iz+1) % Nz;
  int nei = ixy + Nxy * (iz2 + Nz * it);

  real_t bc2 = 1.0;
  if(iz == Nz-1) bc2 = bc_z;

  real_t vt[NVC*ND], vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
      vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
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
      ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
      ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] +=  wt1[1];
    v2[IDX2_SP_I(ic, 2, ist)] += -wt1[0];
    v2[IDX2_SP_R(ic, 3, ist)] += -wt2[1];
    v2[IDX2_SP_I(ic, 3, ist)] +=  wt2[0];
  }

}

//====================================================================
void mult_wilson_zpb(real_t* v2, real_t* u, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int bc_z = bc[2];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  mult_wilson_zpb_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                      Nx, Ny, Nz, Nt, bc_z);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_zm1_kernel(real_t* buf, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt)
{
  int Nxy = Nx * Ny;

  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int ixy = ixyt % Nxy;
  int it  = ixyt / Nxy;

  int iz = Nz-1;
  int ist  = ixy + Nxy * (iz + Nz * it);

  real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
      vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
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
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
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

//====================================================================
void mult_wilson_zm1(real_t* buf, real_t* u, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxyt = Nx * Ny * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nxyt/nth;
  mult_wilson_zm1_kernel<<<nbl,nth>>>(buf_dev, u_dev, v1_dev,
				      Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_zm2_kernel(real_t* v2, real_t* buf, 
                            int Nx, int Ny, int Nz, int Nt, int bc_z)
{
  int Nxy = Nx * Ny;

  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int ixy = ixyt % Nxy;
  int it  = ixyt / Nxy;
  int iz = 0;
  int ist  = ixy + Nxy * (iz + Nz * it);

  real_t bc2 = bc_z;

  real_t wt1[2], wt2[2];

  for(int ic = 0; ic < NC; ++ic){
    int icr = 2 * ic;
    int ici = 2 * ic + 1;

    wt1[0] = bc2 * buf[icr       + 2*NVC*ixyt];
    wt1[1] = bc2 * buf[ici       + 2*NVC*ixyt];
    wt2[0] = bc2 * buf[icr + NVC + 2*NVC*ixyt];
    wt2[1] = bc2 * buf[ici + NVC + 2*NVC*ixyt];

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] += -wt1[1];
    v2[IDX2_SP_I(ic, 2, ist)] +=  wt1[0];
    v2[IDX2_SP_R(ic, 3, ist)] +=  wt2[1];
    v2[IDX2_SP_I(ic, 3, ist)] += -wt2[0];
  }

}

//====================================================================
void mult_wilson_zm2(real_t* v2, real_t* buf,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxyt = Nx * Ny * Nt;

  int bc_z = bc[2];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);

  int nth = VECTOR_LENGTH;
  int nbl = Nxyt/nth;
  mult_wilson_zm2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
				      Nx, Ny, Nz, Nt, bc_z);

  CHECK(cudaDeviceSynchronize());

}

//==================================================================== 
__global__
void mult_wilson_zmb_kernel(real_t* v2, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_z)
{
  int Nxy = Nx * Ny;

  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ixy = ist % Nxy;
  int iz  = (ist/Nxy) % Nz;
  int it  = ist/(Nxy*Nz);
  int iz2 = (iz-1+Nz) % Nz;
  int nei = ixy + Nxy * (iz2 + Nz * it);
  real_t bc2 = 1.0;
  if(iz == 0) bc2 = bc_z;

  real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  for(int id = 0; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
     int icr = 2*ic;
     int ici = 2*ic + 1;
     vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
     vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
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
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, nei)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, nei)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  bc2 * wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  bc2 * wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  bc2 * wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  bc2 * wt2[1];
    v2[IDX2_SP_R(ic, 2, ist)] += -bc2 * wt1[1];
    v2[IDX2_SP_I(ic, 2, ist)] +=  bc2 * wt1[0];
    v2[IDX2_SP_R(ic, 3, ist)] +=  bc2 * wt2[1];
    v2[IDX2_SP_I(ic, 3, ist)] += -bc2 * wt2[0];
  }

}

//====================================================================
void mult_wilson_zmb(real_t* v2, real_t* u, real_t* v1,
                     int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int bc_z = bc[2];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  mult_wilson_zmb_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                      Nx, Ny, Nz, Nt, bc_z);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_tp1_dirac_kernel(real_t* buf, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_t)
{
  int Nxyz = Nx * Ny * Nz;

  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int it  = 0;
  int ist = ixyz + Nxyz * it;

  real_t bc2 = bc_t;

  real_t vt[NVC*ND];

  for(int id = ND2; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
      vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
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

//====================================================================
void mult_wilson_tp1_dirac(real_t* buf, real_t* v1,
                           int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxyz = Nx * Ny * Nz;

  int bc_t = bc[3];

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nxyz/nth;
  mult_wilson_tp1_dirac_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
                                            Nx, Ny, Nz, Nt, bc_t);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_tp2_dirac_kernel(real_t* v2, real_t* u, real_t* buf, 
                                  int Nx, int Ny, int Nz, int Nt)
{
  int Nxyz = Nx * Ny * Nz;

  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

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
      ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
      ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 2, ist)] += wt1[0];
    v2[IDX2_SP_I(ic, 2, ist)] += wt1[1];
    v2[IDX2_SP_R(ic, 3, ist)] += wt2[0];
    v2[IDX2_SP_I(ic, 3, ist)] += wt2[1];
  }

}

//====================================================================
void mult_wilson_tp2_dirac(real_t* v2,  real_t* u, real_t* buf,
                           int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxyz = Nx * Ny * Nz;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* buf_dev = (real_t*)dev_ptr(buf);

  int nth = VECTOR_LENGTH;
  int nbl = Nxyz/nth;
  mult_wilson_tp2_dirac_kernel<<<nbl,nth>>>(v2_dev, u_dev, buf_dev,
                                            Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_tpb_dirac_kernel(real_t* v2, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_t)
{
  int Nxyz = Nx * Ny * Nz;

  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ixyz = ist % Nxyz;
  int it  = ist/Nxyz;

  int it2 = (it+1) % Nt;
  int nei = ixyz + Nxyz * it2;

  real_t bc2 = 1.0;
  if(it == Nt-1) bc2 = bc_t;

  real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  for(int id = ND2; id < ND; ++id){
    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
      vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
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
      ut[2*ic2  ] = u[IDX2_G_R(ic2, ic, ist)];
      ut[2*ic2+1] = u[IDX2_G_I(ic2, ic, ist)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 2, ist)] += wt1[0];
    v2[IDX2_SP_I(ic, 2, ist)] += wt1[1];
    v2[IDX2_SP_R(ic, 3, ist)] += wt2[0];
    v2[IDX2_SP_I(ic, 3, ist)] += wt2[1];
  }

}

//====================================================================
void mult_wilson_tpb_dirac(real_t* v2, real_t* u, real_t* v1,
                           int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int bc_t = bc[3];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  mult_wilson_tpb_dirac_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                            Nx, Ny, Nz, Nt, bc_t);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_tm1_dirac_kernel(real_t* buf, real_t* u, real_t* v1,
                                  int Nx, int Ny, int Nz, int Nt)
{
  int Nxyz = Nx * Ny * Nz;

  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int it   = Nt-1;
  int ist  = ixyz + Nxyz * it;

  real_t vt[2*NC*ND], vt1[NVC], vt2[NVC], ut[NDF];

  for(int id = 0; id < ND2; ++id){
    for(int ic = 0; ic < NC; ++ic){
      vt[2*ic   + NVC*id] = v1[ IDX2_SP_R(ic, id, ist) ];
      vt[2*ic+1 + NVC*id] = v1[ IDX2_SP_I(ic, id, ist) ];
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
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, ist)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, ist)];
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

//====================================================================
void mult_wilson_tm1_dirac(real_t* buf, real_t* u, real_t* v1,
                           int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxyz = Nx * Ny * Nz;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nxyz/nth;
  mult_wilson_tm1_dirac_kernel<<<nbl,nth>>>(buf_dev, u_dev, v1_dev,
                                            Nx, Ny, Nz, Nt);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
__global__
void mult_wilson_tm2_dirac_kernel(real_t* v2, real_t* buf, 
                           int Nx, int Ny, int Nz, int Nt, int bc_t)
{
  int Nxyz = Nx * Ny * Nz;

  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int it  = 0;
  int ist = ixyz + Nxyz * it;

  real_t bc2 = bc_t;

  real_t wt1[2], wt2[2];

  for(int ic = 0; ic < NC; ++ic){
    int icr = 2 * ic;
    int ici = 2 * ic + 1;

    wt1[0] = bc2 * buf[icr       + 2 * NVC * ixyz];
    wt1[1] = bc2 * buf[ici       + 2 * NVC * ixyz];
    wt2[0] = bc2 * buf[icr + NVC + 2 * NVC * ixyz];
    wt2[1] = bc2 * buf[ici + NVC + 2 * NVC * ixyz];

    v2[IDX2_SP_R(ic, 0, ist)] +=  wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  wt2[1];
  }

}

//====================================================================
void mult_wilson_tm2_dirac(real_t* v2, real_t* buf,
                           int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nxyz = Nx * Ny * Nz;

  int bc_t = bc[3];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);

  int nth = VECTOR_LENGTH;
  int nbl = Nxyz/nth;
  mult_wilson_tm2_dirac_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
                                            Nx, Ny, Nz, Nt, bc_t);

  CHECK(cudaDeviceSynchronize());

}

//==================================================================== 
__global__
void mult_wilson_tmb_dirac_kernel(real_t* v2, real_t* u, real_t* v1,
                            int Nx, int Ny, int Nz, int Nt, int bc_t)
{
  int Nxyz = Nx * Ny * Nz;

  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ixyz = ist % Nxyz;
  int it   = ist/Nxyz;

  int it2 = (it-1+Nt) % Nt;
  int nei = ixyz + Nxyz * it2;

  real_t bc2 = 1.0;
  if(it == 0) bc2 = bc_t;

  real_t vt[2*NC*ND];
  real_t vt1[NVC], vt2[NVC], ut[NDF];
  real_t wt1[2], wt2[2];

  for(int id = 0; id < ND2; ++id){
    for(int ic = 0; ic < NC; ++ic){
      int icr = 2*ic;
      int ici = 2*ic + 1;
      vt[icr + NVC*id] = v1[ IDX2_SP_R(ic, id, nei) ];
      vt[ici + NVC*id] = v1[ IDX2_SP_I(ic, id, nei) ];
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
      ut[2*ic2  ] =   u[IDX2_G_R(ic, ic2, nei)];
      ut[2*ic2+1] = - u[IDX2_G_I(ic, ic2, nei)];
    }

    wt1[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt1[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt1[0],vt1[1],vt1[2],vt1[3],vt1[4],vt1[5]);
    wt2[0] = MULT_UV_R(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);
    wt2[1] = MULT_UV_I(ut[0], ut[1], ut[2], ut[3], ut[4], ut[5],
                      vt2[0],vt2[1],vt2[2],vt2[3],vt2[4],vt2[5]);

    v2[IDX2_SP_R(ic, 0, ist)] +=  bc2 * wt1[0];
    v2[IDX2_SP_I(ic, 0, ist)] +=  bc2 * wt1[1];
    v2[IDX2_SP_R(ic, 1, ist)] +=  bc2 * wt2[0];
    v2[IDX2_SP_I(ic, 1, ist)] +=  bc2 * wt2[1];
  }

}

//====================================================================
void mult_wilson_tmb_dirac(real_t* v2, real_t* u, real_t* v1,
                           int *Nsize, int *bc, int Nc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int bc_t = bc[3];

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* u_dev   = (real_t*)dev_ptr(u);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  mult_wilson_tmb_dirac_kernel<<<nbl,nth>>>(v2_dev, u_dev, v1_dev,
                                            Nx, Ny, Nz, Nt, bc_t);

  CHECK(cudaDeviceSynchronize());

}

//====================================================================
void mult_wilson_tp1_chiral(real_t* buf, real_t* v1,
                            int *Nsize, int *bc, int Nc)
{
  printf("mult_wilson_tp1_chiral is not implemented\n");
  exit(1);
}

//====================================================================
void mult_wilson_tm1_chiral(real_t* buf, real_t* u, real_t* v1,
                            int *Nsize, int *bc, int Nc)
{
  printf("mult_wilson_tm1_chiral is not implemented\n");
  exit(1);
}

//====================================================================
void mult_wilson_tp2_chiral(real_t* v2,  real_t* u, real_t* buf,
                            int *Nsize, int *bc, int Nc)
{
  printf("mult_wilson_tp2_chiral is not implemented\n");
  exit(1);
}

//====================================================================
void mult_wilson_tm2_chiral(real_t* v2, real_t* buf,
                            int *Nsize, int *bc, int Nc)
{
  printf("mult_wilson_tm2_chiral is not implemented\n");
  exit(1);
}

//====================================================================
void mult_wilson_tpb_chiral(real_t* v2, real_t* u, real_t* v1,
                            int *Nsize, int *bc, int Nc)
{
  printf("mult_wilson_tpb_chiral is not implemented\n");
  exit(1);
}

//====================================================================
void mult_wilson_tmb_chiral(real_t* v2, real_t* u, real_t* v1,
                            int *Nsize, int *bc, int Nc)
{
  printf("mult_wilson_tmb_chiral is not implemented\n");
  exit(1);
}

//============================================================END=====

