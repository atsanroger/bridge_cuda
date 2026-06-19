/*!

        @file    aligned_allocator_impl.h

        @brief   allocator with alignment

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $

        @version $LastChangedRevision: 2499 $
*/
//====================================================================

#ifndef __ALINGED_ALLOCATOR_IMPL_H__
#define __ALINGED_ALLOCATOR_IMPL_H__

#include <stdlib.h>
#include <cstddef>

/*
 cf.
  Programing Language C++ 4th eidtion
  Bjarne Stroustrup
  sec. 34.4

see also
  allocator.h in stdc++ lib

 */
template<typename _Tp, int AlignmentSize>
struct aligned_allocator_impl
{
  typedef size_t      size_type;
  typedef ptrdiff_t   difference_type;
  typedef _Tp *       pointer;
  typedef const _Tp * const_pointer;
  typedef _Tp&        reference;
  typedef const _Tp&  const_reference;
  typedef _Tp         value_type;

  aligned_allocator_impl() throw() {}
  ~aligned_allocator_impl() throw() {}

  template<typename T>
  using allocator_t = aligned_allocator_impl<T, AlignmentSize>;
  template<typename _Tp1>
  struct rebind
  { typedef allocator_t<_Tp1> other; };

  _Tp *allocate(size_t n)
  {
    void *p       = 0;
    void **memptr = &p;
    int  err      = posix_memalign(memptr, AlignmentSize, sizeof(_Tp) * n);
    if (err) {
      throw std::bad_alloc();
    }
    return reinterpret_cast<_Tp *>(p);
  }

  void deallocate(_Tp *p, size_t)
  {
    free(p);
  }
};

template<typename _T1, int S1, typename _T2, int S2>
inline bool
operator==(const aligned_allocator_impl<_T1, S1>&,
           const aligned_allocator_impl<_T2, S2>&)
{ return true; }

template<typename _Tp, int S>
inline bool
operator==(const aligned_allocator_impl<_Tp, S>&,
           const aligned_allocator_impl<_Tp, S>&)
{ return true; }

template<typename _T1, int S1, typename _T2, int S2>
inline bool
operator!=(const aligned_allocator_impl<_T1, S1>&,
           const aligned_allocator_impl<_T2, S2>&)
{ return false; }

template<typename _Tp, int S>
inline bool
operator!=(const aligned_allocator_impl<_Tp, S>&, const aligned_allocator_impl<_Tp, S>&)
{ return false; }



template<typename _Tp, int AlignmentSize, int OffsetSize>
struct aligned_allocator_offset_impl
{
  typedef size_t     size_type;
  typedef ptrdiff_t  difference_type;
  typedef _Tp *      pointer;
  typedef const _Tp *const_pointer;
  typedef _Tp&       reference;
  typedef const _Tp& const_reference;
  typedef _Tp        value_type;

  aligned_allocator_offset_impl() throw() {}
  ~aligned_allocator_offset_impl() throw() {}

  template<typename T>
  using allocator_t = aligned_allocator_offset_impl<T, AlignmentSize, OffsetSize>;
  template<typename _Tp1>
  struct rebind
  { typedef allocator_t<_Tp1> other; };

  _Tp *allocate(size_t n)
  {
    void *p       = 0;
    void **memptr = &p;
    assert(AlignmentSize % sizeof(_Tp) == 0);
    int offset = OffsetSize / sizeof(_Tp);
    int err    = posix_memalign(memptr, AlignmentSize, sizeof(_Tp) * (n + 2 * offset));
    if (err) {
      throw std::bad_alloc();
    }
    return reinterpret_cast<_Tp *>(p) + offset;
  }

  void deallocate(_Tp *p, size_t)
  {
    int offset = OffsetSize / sizeof(_Tp);
    free(p - offset);
  }
};

template<typename _T1, int S1, int O1, typename _T2, int S2, int O2>
inline bool
operator==(const aligned_allocator_offset_impl<_T1, S1, O1>&,
           const aligned_allocator_offset_impl<_T2, S2, O2>&)
{ return true; }

template<typename _Tp, int S, int O>
inline bool
operator==(const aligned_allocator_offset_impl<_Tp, S, O>&,
           const aligned_allocator_offset_impl<_Tp, S, O>&)
{ return true; }

template<typename _T1, int S1, int O1, typename _T2, int S2, int O2>
inline bool
operator!=(const aligned_allocator_offset_impl<_T1, S1, O1>&,
           const aligned_allocator_offset_impl<_T2, S2, O2>&)
{ return false; }

template<typename _Tp, int S, int O>
inline bool
operator!=(const aligned_allocator_offset_impl<_Tp, S, O>&, const aligned_allocator_offset_impl<_Tp, S, O>&)
{ return false; }

#endif
