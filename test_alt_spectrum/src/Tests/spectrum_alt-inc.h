/*!
      @file    spectrum_alt-inc.h
      @brief
      @author  Hideo Matsufuru
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 12:31:44 #$
      @version $LastChangedRevision: 2500 $
*/

#include "spectrum_alt.h"
#include "lib/Tools/timer.h"

//====================================================================
template<typename AFIELD>
void Spectrum_alt::check_linearAlgebra(int Nin, int Nvol, int Nex,
                                       int Nvec)
{
  int Nrepeat = 500;
  typedef typename AFIELD::real_t real_t;

  vout.general("\n");
  vout.general("Linear algebra performance\n");
  vout.general("  Nin  = %d\n", Nin);
  vout.general("  Nvol = %d\n", Nvol);
  vout.general("  Nex  = %d\n", Nex);
  vout.general("  Nvec = %d\n", Nvec);
  if (sizeof(real_t) == 4) {
    vout.general("  precision: single\n");
  } else if (sizeof(real_t) == 8) {
    vout.general("  precision: double\n");
  }

  vector<AFIELD> x(Nvec), y(Nvec);
  vector<real_t> a(Nvec);
  for (int i = 0; i < Nvec; ++i) {
    x[i].reset(Nin, Nvol, Nex);
    y[i].reset(Nin, Nvol, Nex);
    real_t val = 0.1 * real_t(i + 1);
    a[i] = 1.0 / val;
    x[i].set(val);
    y[i].set(0.0);
  }

  double flop_opr, flop_total, elapsed_time, flops, gflops;

  unique_ptr<Timer> timer(new Timer);

  // axpy
  timer->start();
#pragma omp parallel
  {
    for (int i = 0; i < Nrepeat; ++i) {
      for (int j = 0; j < Nvec; ++j) {
        axpy(y[j], a[j], x[j]);
      }
    }
  }
  timer->stop();

  flop_opr = 2.0 * double(Nin * Nvol * Nex);  // axpy

  flop_total   = flop_opr * double(Nvec * Nrepeat);
  elapsed_time = timer->elapsed_sec();
  flops        = flop_total / elapsed_time;
  gflops       = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Performance of axpy:\n");
  vout.general(m_vl, "  Elapsed time = %12.6f sec\n", elapsed_time);
  vout.general(m_vl, "  Flop(axpy)   = %f\n", flop_total);
  vout.general(m_vl, "  Performance  = %12.3f GFlops\n", gflops);

  timer->reset();

  // dot
  timer->start();

  for (int j = 0; j < Nvec; ++j) {
    a[j] = 0.0;
  }

#pragma omp parallel
  {
    for (int i = 0; i < Nrepeat; ++i) {
      for (int j = 0; j < Nvec; ++j) {
        double a2 = dot(y[j], x[j]);
#pragma omp master
        a[j] += a2;
#pragma omp barrier
      }
    }
  }
  timer->stop();

  // check of result.
  vout.general(m_vl, "\n");
  vout.general(m_vl, "check of result:\n");
  real_t res = 0.0;
  for (int j = 0; j < Nvec; ++j) {
    a[j] = a[j] / (Nrepeat * Nrepeat);
    real_t a_exp = 0.1 * (j + 1) * x[j].size();
    // vout.general(m_vl, "a[%d] = %f (%f) \n", j, a[j], a_exp);
    real_t res1 = a[j] - a_exp;
    res1 = sqrt(res1 * res1);
    res += res1;
  }
  if (res > 1.e-7)
    vout.general(m_vl, "worning: too large difference!\n", res);

  flop_opr = 2.0 * double(Nin * Nvol * Nex);

  flop_total   = flop_opr * double(Nvec * Nrepeat);
  elapsed_time = timer->elapsed_sec();
  flops        = flop_total / elapsed_time;
  gflops       = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Performance of dot:\n");
  vout.general(m_vl, "  Elapsed time = %12.6f sec\n", elapsed_time);
  vout.general(m_vl, "  Flop(axpy)   = %f\n", flop_total);
  vout.general(m_vl, "  Performance  = %12.3f GFlops\n", gflops);

  timer->reset();

  // norm2
#pragma omp parallel
  {
    timer->start();
    for (int i = 0; i < Nrepeat; ++i) {
      for (int j = 0; j < Nvec; ++j) {
        double a = y[j].norm2();
      }
    }
  }
  timer->stop();

  flop_opr = 2.0 * double(Nin * Nvol * Nex);

  flop_total   = flop_opr * double(Nvec * Nrepeat);
  elapsed_time = timer->elapsed_sec();
  flops        = flop_total / elapsed_time;
  gflops       = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "Performance of norm2:\n");
  vout.general(m_vl, "  Elapsed time = %12.6f sec\n", elapsed_time);
  vout.general(m_vl, "  Flop(axpy)   = %f\n", flop_total);
  vout.general(m_vl, "  Performance  = %12.3f GFlops\n", gflops);
}


//============================================================END=====
