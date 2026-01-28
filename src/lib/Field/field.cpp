/*!
        @file    field.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Field/field.h"

#include <cstring>
using std::string;

#include "ResourceManager/threadManager.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

#include "Field/field_thread-inc.h"

const std::string Field::class_name = "Field";

//====================================================================
void Field::check()
{
  // ThreadManager::assert_single_thread(class_name);
  //    vout.general("Field was constructed.\n");
}


//====================================================================
void Field::set(double a)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  double *yp = this->ptr(0);

  for (int ex = 0; ex < m_Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int kv = m_Nin * (site + m_Nvol * ex);
      for (int in = 0; in < m_Nin; ++in) {
        yp[in + kv] = a;
      }
    }
  }
}


//====================================================================
void Field::setc(const double a)
{
  if (m_element_type == Element_type::REAL) {
    set(a);
  } else if (m_element_type == Element_type::COMPLEX) {
    dcomplex c = cmplx(a, 0.0);
    setc(c);
  } else {
    vout.crucial("Error at %s: unsupported arg types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Field::setc(const dcomplex a)
{
  if (imag(a) == 0.0) {
    return set(real(a));
  } else if (m_element_type == Element_type::COMPLEX) {
    double *yp  = this->ptr(0);
    double ar   = real(a);
    double ai   = imag(a);
    int    Nin2 = m_Nin / 2;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    for (int ex = 0; ex < m_Nex; ++ex) {
      for (int site = is; site < ns; ++site) {
        int kv = m_Nin * (site + m_Nvol * ex);
        for (int k = 0; k < Nin2; ++k) {
          int kr = 2 * k;
          int ki = 2 * k + 1;
          yp[kr + kv] = ar;
          yp[ki + kv] = ai;
        }
      }
    }
  } else if (m_element_type == Element_type::REAL) {
    double ar = real(a);
    double ai = imag(a);

    if (fabs(ai) < fabs(ar) * CommonParameters::epsilon_criterion()) {
      return set(real(a));
    } else {
      vout.crucial("Error at %s: real vector and complex parameter.\n",
                   __func__);
      exit(EXIT_FAILURE);
    }
  } else {
    vout.crucial("Error at %s: unsupported arg types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
double Field::norm2() const
{
  const double *yp = this->ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  double a = 0.0;

  for (int ex = 0; ex < m_Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int kv = m_Nin * (site + m_Nvol * ex);
      for (int in = 0; in < m_Nin; ++in) {
        a += yp[in + kv] * yp[in + kv];
      }
    }
  }

  ThreadManager::reduce_sum_global(a, ith, nth);

  return a;
}


//====================================================================
void Field::xI()
{
  if (m_element_type == Element_type::REAL) {
    vout.crucial("Error at %s: xI() is not available for real field.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  double *yp  = this->ptr(0);
  int    Nin2 = m_Nin / 2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  for (int ex = 0; ex < m_Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int kv = m_Nin * (site + m_Nvol * ex);
      for (int k = 0; k < Nin2; ++k) {
        int    kr = 2 * k;
        int    ki = 2 * k + 1;
        double yr = yp[kr + kv];
        double yi = yp[ki + kv];
        yp[kr + kv] = -yi;
        yp[ki + kv] = yr;
      }
    }
  }
}


//====================================================================
void Field::stat(double& Fave, double& Fmax, double& Fdev) const
{
#pragma omp barrier

  double sum  = 0.0;
  double sum2 = 0.0;
  double vmax = 0.0;

  const double *yp = this->ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  vmax = 0.0;
  for (int ex = 0; ex < m_Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      double fst = 0.0;
      for (int in = 0; in < m_Nin; ++in) {
        double fv = yp[myindex(in, site, ex)];
        fst += fv * fv;
      }
      sum2 += fst;
      fst   = sqrt(fst);
      sum  += fst;
      if (fst > vmax) vmax = fst;
    }
  }

  ThreadManager::reduce_sum_global(sum, ith, nth);
  ThreadManager::reduce_sum_global(sum2, ith, nth);
  ThreadManager::reduce_max_global(vmax, ith, nth);

  double vfac = double(m_Nvol) * double(m_Nex)
                * double(CommonParameters::NPE());
  Fave = sum / vfac;
  Fdev = sqrt(sum2 / vfac - Fave * Fave);
  Fmax = vmax;

#pragma omp barrier
}


//====================================================================
void copy(Field& y, const Field& x)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  int Nex  = y.nex();
  assert(x.check_size(Nin, Nvol, Nex));

  double       *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int kv = Nin * (site + Nvol * ex);
      for (int in = 0; in < Nin; ++in) {
        yp[in + kv] = xp[in + kv];
      }
    }
  }
}


//====================================================================
void copy(Field& y, const int exy, const Field& x, const int exx)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();

  assert(x.nin() == Nin);
  assert(x.nvol() == Nvol);

  double       *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int site = is; site < ns; ++site) {
    int kv = Nin * site;
    for (int in = 0; in < Nin; ++in) {
      yp[in + kv] = xp[in + kv];
    }
  }
}


//====================================================================
void scal(Field& x, const double a)
{
  int Nin  = x.nin();
  int Nvol = x.nvol();
  int Nex  = x.nex();

  double *xp = x.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int kv = Nin * (site + Nvol * ex);
      for (int in = 0; in < Nin; ++in) {
        xp[in + kv] *= a;
      }
    }
  }
}


//====================================================================
void scal(Field& x, const int exx, const double a)
{
  int Nin  = x.nin();
  int Nvol = x.nvol();

  double *xp = x.ptr(0, 0, exx);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int site = is; site < ns; ++site) {
    int kv = Nin * site;
    for (int in = 0; in < Nin; ++in) {
      xp[in + kv] *= a;
    }
  }
}


//====================================================================
void scal(Field& x, const dcomplex a)
{
  if (imag(a) == 0.0) {
    return scal(x, real(a));
  } else if (x.field_element_type() == Element_type::COMPLEX) {
    int Nin  = x.nin();
    int Nvol = x.nvol();
    int Nex  = x.nex();

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nvol);

    double *xp  = x.ptr(0);
    double ar   = real(a);
    double ai   = imag(a);
    int    Nin2 = Nin / 2;

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = is; site < ns; ++site) {
        int kv = Nin * (site + Nvol * ex);
        for (int k = 0; k < Nin2; ++k) {
          int    kr = 2 * k;
          int    ki = 2 * k + 1;
          double xr = xp[kr + kv];
          double xi = xp[ki + kv];
          xp[kr + kv] = ar * xr - ai * xi;
          xp[ki + kv] = ar * xi + ai * xr;
        }
      }
    }
  } else {
    vout.crucial("Error at %s: real vector and complex parameter.\n",
                 __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void scal(Field& x, const int exx, const dcomplex a)
{
  if (imag(a) == 0.0) {
    return scal(x, exx, real(a));
  } else if (x.field_element_type() == Element_type::COMPLEX) {
    int Nin  = x.nin();
    int Nvol = x.nvol();
    int Nex  = x.nex();

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nvol);

    double *xp  = x.ptr(0, 0, exx);
    double ar   = real(a);
    double ai   = imag(a);
    int    Nin2 = Nin / 2;

    for (int site = is; site < ns; ++site) {
      int kv = Nin * site;
      for (int k = 0; k < Nin2; ++k) {
        int    kr = 2 * k;
        int    ki = 2 * k + 1;
        double xr = xp[kr + kv];
        double xi = xp[ki + kv];
        xp[kr + kv] = ar * xr - ai * xi;
        xp[ki + kv] = ar * xi + ai * xr;
      }
    }
  } else {
    vout.crucial("Error at %s: real vector and complex parameter.\n",
                 __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void axpy(Field& y, const double a, const Field& x)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  int Nex  = y.nex();
  assert(x.check_size(Nin, Nvol, Nex));

  double       *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int kv = Nin * (site + Nvol * ex);
      for (int in = 0; in < Nin; ++in) {
        yp[in + kv] += a * xp[in + kv];
      }
    }
  }
}


//====================================================================
void axpy(Field& y, const int exy,
          const double a, const Field& x, const int exx)
{
  assert(x.nin() == y.nin());
  assert(x.nvol() == y.nvol());

  int Nin  = y.nin();
  int Nvol = y.nvol();

  double       *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int site = is; site < ns; ++site) {
    int kv = Nin * site;
    for (int in = 0; in < Nin; ++in) {
      yp[in + kv] += a * xp[in + kv];
    }
  }
}


//====================================================================
void axpy(Field& y, const dcomplex a, const Field& x)
{
  if (imag(a) == 0.0) {
    return axpy(y, real(a), x);
  } else if ((y.field_element_type() == Element_type::COMPLEX) &&
             (x.field_element_type() == Element_type::COMPLEX)) {
    int Nin  = y.nin();
    int Nvol = y.nvol();
    int Nex  = y.nex();
    assert(x.check_size(Nin, Nvol, Nex));

    double       *yp = y.ptr(0);
    const double *xp = x.ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nvol);

    double ar   = real(a);
    double ai   = imag(a);
    int    Nin2 = Nin / 2;

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = is; site < ns; ++site) {
        int kv = Nin * (site + Nvol * ex);
        for (int k = 0; k < Nin2; ++k) {
          int kr = 2 * k;
          int ki = 2 * k + 1;
          yp[kr + kv] += ar * xp[kr + kv] - ai * xp[ki + kv];
          yp[ki + kv] += ar * xp[ki + kv] + ai * xp[kr + kv];
        }
      }
    }
  } else {
    vout.crucial("Error at %s: unsupported types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void axpy(Field& y, const int exy,
          const dcomplex a, const Field& x, const int exx)
{
  if (imag(a) == 0.0) {
    return axpy(y, exy, real(a), x, exx);
  } else if ((y.field_element_type() == Element_type::COMPLEX) &&
             (x.field_element_type() == Element_type::COMPLEX)) {
    int Nin  = y.nin();
    int Nvol = y.nvol();
    assert(x.nin() == Nin);
    assert(x.nvol() == Nvol);

    double       *yp = y.ptr(0, 0, exy);
    const double *xp = x.ptr(0, 0, exx);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nvol);

    double ar   = real(a);
    double ai   = imag(a);
    int    Nin2 = Nin / 2;

    for (int site = is; site < ns; ++site) {
      int kv = Nin * site;
      for (int k = 0; k < Nin2; ++k) {
        int kr = 2 * k;
        int ki = 2 * k + 1;
        yp[kr + kv] += ar * xp[kr + kv] - ai * xp[ki + kv];
        yp[ki + kv] += ar * xp[ki + kv] + ai * xp[kr + kv];
      }
    }
  } else {
    vout.crucial("Error at %s: unsupported types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void aypx(const double a, Field& y, const Field& x)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  int Nex  = y.nex();
  assert(x.check_size(Nin, Nvol, Nex));

  double       *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int kv = Nin * (site + Nvol * ex);
      for (int in = 0; in < Nin; ++in) {
        yp[in + kv] = a * yp[in + kv] + xp[in + kv];
      }
    }
  }
}


//====================================================================
void aypx(const dcomplex a, Field& y, const Field& x)
{
  if (imag(a) == 0.0) {
    return aypx(real(a), y, x);
  } else if ((y.field_element_type() == Element_type::COMPLEX) &&
             (x.field_element_type() == Element_type::COMPLEX)) {
    int Nin  = y.nin();
    int Nvol = y.nvol();
    int Nex  = y.nex();
    assert(x.check_size(Nin, Nvol, Nex));

    double       *yp = y.ptr(0);
    const double *xp = x.ptr(0);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nvol);

    double ar   = real(a);
    double ai   = imag(a);
    int    Nin2 = Nin / 2;

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = is; site < ns; ++site) {
        int kv = Nin * (site + Nvol * ex);
        for (int k = 0; k < Nin2; ++k) {
          int    kr = 2 * k;
          int    ki = 2 * k + 1;
          double yr = yp[kr + kv];
          double yi = yp[ki + kv];
          yp[kr + kv] = ar * yr - ai * yi + xp[kr + kv];
          yp[ki + kv] = ar * yi + ai * yr + xp[ki + kv];
        }
      }
    }
  } else {
    vout.crucial("Error at %s: unsupported types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
double dot(const Field& y, const Field& x)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  int Nex  = y.nex();
  assert(x.check_size(Nin, Nvol, Nex));

  const double *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double a = 0.0;

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int kv = Nin * (site + Nvol * ex);
      for (int in = 0; in < Nin; ++in) {
        a += yp[in + kv] * xp[in + kv];
      }
    }
  }

  ThreadManager::reduce_sum_global(a, ith, nth);

  return a;
}


//====================================================================
double dot(const Field& y, const int exy, const Field& x, const int exx)
{
  assert(x.nin() == y.nin());
  assert(x.nvol() == y.nvol());

  int Nin  = y.nin();
  int Nvol = y.nvol();

  const double *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double a = 0.0;

  for (int site = is; site < ns; ++site) {
    int kv = Nin * site;
    for (int in = 0; in < Nin; ++in) {
      a += yp[in + kv] * xp[in + kv];
    }
  }

  ThreadManager::reduce_sum_global(a, ith, nth);

  return a;
}


//====================================================================
void dot_and_norm2(double& yx, double& y2, double& x2,
                   const Field& y, const int exy,
                   const Field& x, const int exx)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  assert(x.nin() == Nin);
  assert(x.nvol() == Nvol);

  const double *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double sum_yx = 0.0;
  double sum_x2 = 0.0;
  double sum_y2 = 0.0;

  for (int site = is; site < ns; ++site) {
    int kv = Nin * site;
    for (int in = 0; in < Nin; ++in) {
      sum_yx += yp[in + kv] * xp[in + kv];
      sum_x2 += xp[in + kv] * xp[in + kv];
      sum_y2 += yp[in + kv] * yp[in + kv];
    }
  }

  double prd[3] = { sum_yx, sum_x2, sum_y2 };
  ThreadManager::reduce_sum_global(prd, 3, ith, nth);
  yx = prd[0];
  y2 = prd[1];
  x2 = prd[2];
}


//====================================================================
void dot_and_norm2(double& yx, double& y2, double& x2,
                   const Field& y, const Field& x)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  int Nex  = y.nex();
  assert(x.check_size(Nin, Nvol, Nex));

  const double *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double sum_yx = 0.0;
  double sum_x2 = 0.0;
  double sum_y2 = 0.0;

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      int kv = Nin * (site + Nvol * ex);
      for (int in = 0; in < Nin; ++in) {
        sum_yx += yp[in + kv] * xp[in + kv];
        sum_x2 += xp[in + kv] * xp[in + kv];
        sum_y2 += yp[in + kv] * yp[in + kv];
      }
    }
  }

  double prd[3] = { sum_yx, sum_x2, sum_y2 };
  ThreadManager::reduce_sum_global(prd, 3, ith, nth);
  yx = prd[0];
  y2 = prd[1];
  x2 = prd[2];
}


//====================================================================
dcomplex dotc(const Field& y, const Field& x)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  int Nex  = y.nex();
  assert(x.check_size(Nin, Nvol, Nex));

  const double *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  if ((y.field_element_type() == Element_type::COMPLEX) &&
      (x.field_element_type() == Element_type::COMPLEX)) {
    double prdr = 0.0;
    double prdi = 0.0;
    int    Nin2 = Nin / 2;

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = is; site < ns; ++site) {
        int kv = Nin * (site + Nvol * ex);
        for (int k = 0; k < Nin2; ++k) {
          int kr = 2 * k;
          int ki = 2 * k + 1;
          prdr += yp[kr + kv] * xp[kr + kv] + yp[ki + kv] * xp[ki + kv];
          prdi += yp[kr + kv] * xp[ki + kv] - yp[ki + kv] * xp[kr + kv];
        }
      }
    }

    double prd[2] = { prdr, prdi };
    ThreadManager::reduce_sum_global(prd, 2, ith, nth);

    return cmplx(prd[0], prd[1]);
  } else if ((y.field_element_type() == Element_type::REAL) &&
             (y.field_element_type() == Element_type::REAL)) {
    return cmplx(dot(y, x), 0.0);
  } else {
    vout.crucial("Error at %s: unsupported arg types\n", __func__);
    exit(EXIT_FAILURE);

    return cmplx(0.0, 0.0);  // never reached.
  }
}


//====================================================================
dcomplex dotc(const Field& y, const int exy, const Field& x, const int exx)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  assert(x.nin() == Nin);
  assert(x.nvol() == Nvol);

  const double *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  if ((y.field_element_type() == Element_type::COMPLEX) &&
      (x.field_element_type() == Element_type::COMPLEX)) {
    double prdr = 0.0;
    double prdi = 0.0;
    int    Nin2 = Nin / 2;

    for (int site = is; site < ns; ++site) {
      int kv = Nin * site;
      for (int k = 0; k < Nin2; ++k) {
        int kr = 2 * k;
        int ki = 2 * k + 1;
        prdr += yp[kr + kv] * xp[kr + kv] + yp[ki + kv] * xp[ki + kv];
        prdi += yp[kr + kv] * xp[ki + kv] - yp[ki + kv] * xp[kr + kv];
      }
    }

    double prd[2] = { prdr, prdi };
    ThreadManager::reduce_sum_global(prd, 2, ith, nth);

    return cmplx(prd[0], prd[1]);
  } else if ((y.field_element_type() == Element_type::REAL) &&
             (y.field_element_type() == Element_type::REAL)) {
    return cmplx(dot(y, exy, x, exx), 0.0);
  } else {
    vout.crucial("Error at %s: unsupported arg types\n", __func__);
    exit(EXIT_FAILURE);

    return cmplx(0.0, 0.0);  // never reached.
  }
}


//====================================================================
void dotc_and_norm2(dcomplex& yx, double& y2, double& x2,
                    const Field& y, const Field& x)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  int Nex  = y.nex();
  assert(x.check_size(Nin, Nvol, Nex));

  const double *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  if ((y.field_element_type() == Element_type::COMPLEX) &&
      (x.field_element_type() == Element_type::COMPLEX)) {
    double prd_r  = 0.0;
    double prd_i  = 0.0;
    double prd_x2 = 0.0;
    double prd_y2 = 0.0;
    int    Nin2   = Nin / 2;

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = is; site < ns; ++site) {
        int kv = Nin * (site + Nvol * ex);
        for (int k = 0; k < Nin2; ++k) {
          int kr = 2 * k;
          int ki = 2 * k + 1;
          prd_r  += yp[kr + kv] * xp[kr + kv] + yp[ki + kv] * xp[ki + kv];
          prd_i  += yp[kr + kv] * xp[ki + kv] - yp[ki + kv] * xp[kr + kv];
          prd_x2 += xp[kr + kv] * xp[kr + kv] + xp[ki + kv] * xp[ki + kv];
          prd_y2 += yp[kr + kv] * yp[kr + kv] + yp[ki + kv] * yp[ki + kv];
        }
      }
    }

    double prd[4] = { prd_r, prd_i, prd_x2, prd_y2 };
    ThreadManager::reduce_sum_global(prd, 4, ith, nth);

    yx = cmplx(prd[0], prd[1]);
    y2 = prd[2];
    x2 = prd[3];
  } else if ((y.field_element_type() == Element_type::REAL) &&
             (y.field_element_type() == Element_type::REAL)) {
    double yx_re = 0.0;
    dot_and_norm2(yx_re, y2, x2, y, x);
    yx = cmplx(yx_re, 0.0);
  } else {
    vout.crucial("Error at %s: unsupported arg types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void dotc_and_norm2(dcomplex& yx, double& y2, double& x2,
                    const Field& y, const int exy,
                    const Field& x, const int exx)
{
  int Nin  = y.nin();
  int Nvol = y.nvol();
  assert(x.nin() == Nin);
  assert(x.nvol() == Nvol);

  const double *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  if ((y.field_element_type() == Element_type::COMPLEX) &&
      (x.field_element_type() == Element_type::COMPLEX)) {
    double prd_r  = 0.0;
    double prd_i  = 0.0;
    double prd_x2 = 0.0;
    double prd_y2 = 0.0;
    int    Nin2   = Nin / 2;

    for (int site = is; site < ns; ++site) {
      int kv = Nin * site;
      for (int k = 0; k < Nin2; ++k) {
        int kr = 2 * k;
        int ki = 2 * k + 1;
        prd_r  += yp[kr + kv] * xp[kr + kv] + yp[ki + kv] * xp[ki + kv];
        prd_i  += yp[kr + kv] * xp[ki + kv] - yp[ki + kv] * xp[kr + kv];
        prd_x2 += xp[kr + kv] * xp[kr + kv] + xp[ki + kv] * xp[ki + kv];
        prd_y2 += yp[kr + kv] * yp[kr + kv] + yp[ki + kv] * yp[ki + kv];
      }
    }

    double prd[4] = { prd_r, prd_i, prd_x2, prd_y2 };
    ThreadManager::reduce_sum_global(prd, 4, ith, nth);

    yx = cmplx(prd[0], prd[1]);
    y2 = prd[2];
    x2 = prd[3];
  } else if ((y.field_element_type() == Element_type::REAL) &&
             (y.field_element_type() == Element_type::REAL)) {
    double yx_re = 0.0;
    dot_and_norm2(yx_re, y2, x2, y, exy, x, exx);
    yx = cmplx(yx_re, 0.0);
  } else {
    vout.crucial("Error at %s: unsupported arg types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//============================================================END=====
