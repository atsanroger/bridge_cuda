/*!
        @file    sorter.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SORTER_INCLUDED
#define SORTER_INCLUDED

#include <cmath>
#include <algorithm>
#include <cassert>
#include <vector>
#include <string>

#include "bridge_complex.h"

//! @class Sorter
//    sort utility class.
//
//  usage:
//    Sorter sorter(sort_type);
//      sort_type =
//        abs_ascending  : absolute ascending order   |v_0| <= |v_1| <= ...
//        abs_descending : absolute descending order  |v_0| >= |v_1| >= ...
//        ascending      : arithmetic ascending order  v_0 <= v_1 <= ...
//        descending     : arithmetic descending order v_0 >= v_1 >= ...
//    sorter.sort(vec);  to sort vec : vector<double>. destructive.
//    sorter.sort_index(vec);  returns index of vec in sorted order.
//    sorter.comp(a, b);  compare two values a, b with the sort condition.
//
// 2022.02.15: I.Kanamori
//   added fcomplex

template<typename T>
class Sorter
{
 public:

  Sorter(const std::string& type);

  virtual ~Sorter();

  //! sort an array of values; v is sorted on exit.
  void sort(std::vector<T>& v);

  //! sort first n elements of an array of values; v is sorted on exit.
  void sort(std::vector<T>& v, const size_t nelem);

  //! sort an array and return list of index; v is sorted on exit.
  std::vector<int> sort_index(std::vector<T>& v);

  //! sort an array and return list of index. first n elemnts are considered.
  std::vector<int> sort_index(std::vector<T>& v, const size_t nelem);

  //! call sort condition.
  bool comp(const T& lhs, const T& rhs);

 private:

  // sort order
  struct by_abs_ascend;
  struct by_abs_descend;
  struct by_ascend;
  struct by_descend;

  struct by_order;
  struct proxy;

  typedef std::pair<int, T> pair_t;

  by_order *m_order;

  // non-copyable
  Sorter(const Sorter&);
  Sorter& operator=(const Sorter&);
};

//====================================================================
//! base class for sort ordering
template<typename T>
struct Sorter<T>::by_order
{
  virtual ~by_order() {}

  virtual
  bool operator()(const T& lhs, const T& rhs) const = 0;

  bool operator()(const pair_t& lhs, const pair_t& rhs) const
  {
    return operator()(lhs.second, rhs.second);
  }
};

//====================================================================
//! absolute ascending order
//  i.e. |v_0| <= |v_1| <= ... <= |v_{n-1}|
template<typename T>
struct Sorter<T>::by_abs_ascend : public Sorter::by_order
{
  bool operator()(const T& lhs, const T& rhs) const
  { return abs(lhs) < abs(rhs); }
};

//====================================================================
//! absolute descending order
//  i.e. |v_0| >= |v_1| >= ... >= |v_{n-1}|
template<typename T>
struct Sorter<T>::by_abs_descend : public Sorter::by_order
{
  bool operator()(const T& lhs, const T& rhs) const
  {
    return abs(lhs) > abs(rhs);
  }
};

//====================================================================
//! ascending order
//  i.e. v_0 <= v_1 <= ... <= v_{n-1}
template<typename T>
struct Sorter<T>::by_ascend : public Sorter::by_order
{
  virtual
  bool operator()(const T& lhs, const T& rhs) const
  {
    return lhs < rhs;
  }
};

//====================================================================
//! descending order
//  i.e. v_0 >= v_1 >= ... >= v_{n-1}
template<typename T>
struct Sorter<T>::by_descend : public Sorter::by_order
{
  virtual
  bool operator()(const T& lhs, const T& rhs) const
  {
    return rhs < lhs;
  }
};

//====================================================================
//! proxy object to pass to stl sort algorithm
template<typename T>
struct Sorter<T>::proxy
{
 private:
  const Sorter::by_order& m_order;

 public:
  proxy(const Sorter::by_order& order) : m_order(order) {}

  bool operator()(const T& lhs, const T& rhs) const
  {
    return m_order(lhs, rhs);
  }

  bool operator()(const pair_t& lhs, const pair_t& rhs) const
  {
    return m_order(lhs, rhs);
  }
};

//====================================================================
//! constructor with sort ordering as a string arg
template<typename T>
Sorter<T>::Sorter(const std::string& type) : m_order(0)
{
  if (type == "abs_ascending") m_order = new by_abs_ascend;
  if (type == "abs_descending") m_order = new by_abs_descend;
  if (type == "ascending") m_order = new by_ascend;
  if (type == "descending") m_order = new by_descend;
}


//====================================================================
//! destructor
template<typename T>
Sorter<T>::~Sorter()
{
  if (m_order) delete m_order;
}


//====================================================================
//! sort an array of values; v is sorted on exit.
template<typename T>
void Sorter<T>::sort(std::vector<T>& v)
{
  assert(m_order != NULL);

  return std::sort(v.begin(), v.end(), proxy(*m_order));
}


//====================================================================
//! sort first n elments of an array of values; v is sorted on exit.
template<typename T>
void Sorter<T>::sort(std::vector<T>& v, const size_t nelem)
{
  assert(m_order != NULL);

  return std::sort(v.begin(), v.begin() + nelem, proxy(*m_order));
}


//====================================================================
//! sort an array and return list of index.
template<typename T>
std::vector<int> Sorter<T>::sort_index(std::vector<T>& v)
{
  assert(m_order != NULL);

  std::vector<pair_t> w(v.size());

  for (size_t i = 0, n = v.size(); i < n; ++i) {
    w[i].first  = i;
    w[i].second = v[i];
  }

  std::sort(w.begin(), w.end(), proxy(*m_order));

  std::vector<int> idx(v.size());

  for (size_t i = 0, n = v.size(); i < n; ++i) {
    idx[i] = w[i].first;
    v[i]   = w[i].second;
  }

  return idx;
}


//====================================================================
//! sort an array and return list of index. first n elements are condidered.
template<typename T>
std::vector<int> Sorter<T>::sort_index(std::vector<T>& v, const size_t nelem)
{
  assert(m_order != NULL);

  std::vector<pair_t> w(v.size());

  for (size_t i = 0, n = v.size(); i < n; ++i) {
    w[i].first  = i;
    w[i].second = v[i];
  }

  std::sort(w.begin(), w.begin() + nelem, proxy(*m_order));

  std::vector<int> idx(v.size());

  for (size_t i = 0, n = v.size(); i < n; ++i) {
    idx[i] = w[i].first;
    v[i]   = w[i].second;
  }

  return idx;
}


//====================================================================
//! call sort condition.
template<typename T>
bool Sorter<T>::comp(const T& lhs, const T& rhs)
{
  assert(m_order != NULL);
  return m_order->operator()(lhs, rhs);
}


//====================================================================
//! ascending order
//  i.e. v_0 <= v_1 <= ... <= v_{n-1}
template<>
struct Sorter<dcomplex>::by_ascend : public Sorter::by_order
{
  virtual
  bool operator()(const dcomplex& lhs, const dcomplex& rhs) const
  {
    vout.crucial("by_ascend: unsupported ordering\n");
    exit(EXIT_FAILURE);
    return false;
  }
};

//====================================================================
//! descending order
//  i.e. v_0 >= v_1 >= ... >= v_{n-1}
template<>
struct Sorter<dcomplex>::by_descend : public Sorter::by_order
{
  virtual
  bool operator()(const dcomplex& lhs, const dcomplex& rhs) const
  {
    vout.crucial("by_ascend: unsupported ordering\n");
    exit(EXIT_FAILURE);
    return false;
  }
};


//====================================================================
//! ascending order
//  i.e. v_0 <= v_1 <= ... <= v_{n-1}
template<>
struct Sorter<fcomplex>::by_ascend : public Sorter::by_order
{
  virtual
  bool operator()(const fcomplex& lhs, const fcomplex& rhs) const
  {
    vout.crucial("by_ascend: unsupported ordering\n");
    exit(EXIT_FAILURE);
    return false;
  }
};

//====================================================================
//! descending order
//  i.e. v_0 >= v_1 >= ... >= v_{n-1}
template<>
struct Sorter<fcomplex>::by_descend : public Sorter::by_order
{
  virtual
  bool operator()(const fcomplex& lhs, const fcomplex& rhs) const
  {
    vout.crucial("by_ascend: unsupported ordering\n");
    exit(EXIT_FAILURE);
    return false;
  }
};


//==========================================================
#endif /* SORTER_INCLUDED */
//==================================================END=====
