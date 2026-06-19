/*!
        @file    pyb_wrap.h

        @brief   pybridge internal template classes

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-02-25 17:26:27 #$

        @version $LastChangedRevision: 2175 $
*/

#ifndef PYWRAP_H_INCLUDED
#define PYWRAP_H_INCLUDED

template <typename T>
class PyHandle
{
 public:
  PyHandle() : ptr_(nullptr) {}
  PyHandle(T* ptr) : ptr_(ptr) {}
  PyHandle(const PyHandle& rhs) : ptr_(rhs.ptr_) {}
  T& operator*() const { return *ptr_; }
  T* operator->() const { return ptr_; }
  T* get() const { return ptr_; }
  void destroy() { delete ptr_; }
  T& operator[](std::size_t idx) const { return ptr_[idx]; }
 private:
  T* ptr_;
};


template <typename T, typename W1>
class
__attribute__((visibility("hidden")))
PyWrap1 : public T
{
 public:
  PyWrap1(py::object obj1)
    : T(obj1.cast<W1*>())
    , obj1_(obj1)
  { ENTER;
    obj1_.inc_ref();
  }
  ~PyWrap1()
  { ENTER;
    obj1_.dec_ref();
  }
 private:
  py::object obj1_;
};

template <typename T, typename W1, typename W2>
class
__attribute__((visibility("hidden")))
PyWrap2 : public T
{
 public:
  PyWrap2(py::object obj1,
          py::object obj2)
    : T(obj1.cast<W1*>(),
        obj2.cast<W2*>())
    , obj1_(obj1)
    , obj2_(obj2)
  { ENTER;
    obj1_.inc_ref();
    obj2_.inc_ref();
  }
  ~PyWrap2()
  { ENTER;
    obj1_.dec_ref();
    obj2_.dec_ref();
  }
 private:
  py::object obj1_;
  py::object obj2_;
};

template <typename T, typename W1, typename W2, typename W3>
class
__attribute__((visibility("hidden")))
PyWrap3 : public T
{
 public:
  PyWrap3(py::object obj1,
          py::object obj2,
          py::object obj3)
    : T(obj1.cast<W1*>(),
        obj2.cast<W2*>(),
        obj3.cast<W3*>())
    , obj1_(obj1)
    , obj2_(obj2)
    , obj3_(obj3)
  { ENTER;
    obj1_.inc_ref();
    obj2_.inc_ref();
    obj3_.inc_ref();
  }
  ~PyWrap3()
  { ENTER;
    obj1_.dec_ref();
    obj2_.dec_ref();
    obj3_.dec_ref();
  }
 private:
  py::object obj1_;
  py::object obj2_;
  py::object obj3_;
};

template <typename T, typename W1, typename W2, typename W3, typename W4>
class
__attribute__((visibility("hidden")))
PyWrap4 : public T
{
 public:
  PyWrap4(py::object obj1,
          py::object obj2,
          py::object obj3,
          py::object obj4)
    : T(obj1.cast<W1*>(),
        obj2.cast<W2*>(),
        obj3.cast<W3*>(),
        obj4.cast<W4*>())
    , obj1_(obj1)
    , obj2_(obj2)
    , obj3_(obj3)
    , obj4_(obj4)
  { ENTER;
    obj1_.inc_ref();
    obj2_.inc_ref();
    obj3_.inc_ref();
    obj4_.inc_ref();
  }
  ~PyWrap4()
  { ENTER;
    obj1_.dec_ref();
    obj2_.dec_ref();
    obj3_.dec_ref();
    obj4_.dec_ref();
  }
 private:
  py::object obj1_;
  py::object obj2_;
  py::object obj3_;
  py::object obj4_;
};

template <typename T, typename W1, typename W2, typename W3, typename W4, typename W5>
class
__attribute__((visibility("hidden")))
PyWrap5 : public T
{
 public:
  PyWrap5(py::object obj1,
          py::object obj2,
          py::object obj3,
          py::object obj4,
          py::object obj5)
    : T(obj1.cast<W1*>(),
        obj2.cast<W2*>(),
        obj3.cast<W3*>(),
        obj4.cast<W4*>(),
        obj5.cast<W5*>())
    , obj1_(obj1)
    , obj2_(obj2)
    , obj3_(obj3)
    , obj4_(obj4)
    , obj5_(obj5)
  { ENTER;
    obj1_.inc_ref();
    obj2_.inc_ref();
    obj3_.inc_ref();
    obj4_.inc_ref();
    obj5_.inc_ref();
  }
  ~PyWrap5()
  { ENTER;
    obj1_.dec_ref();
    obj2_.dec_ref();
    obj3_.dec_ref();
    obj4_.dec_ref();
    obj5_.dec_ref();
  }
 private:
  py::object obj1_;
  py::object obj2_;
  py::object obj3_;
  py::object obj4_;
  py::object obj5_;
};

template <typename T, typename W1, typename W2, typename W3, typename W4, typename W5, typename W6>
class
__attribute__((visibility("hidden")))
PyWrap6 : public T
{
 public:
  PyWrap6(py::object obj1,
          py::object obj2,
          py::object obj3,
          py::object obj4,
          py::object obj5,
          py::object obj6)
    : T(obj1.cast<W1*>(),
        obj2.cast<W2*>(),
        obj3.cast<W3*>(),
        obj4.cast<W4*>(),
        obj5.cast<W5*>(),
        obj6.cast<W6*>())
    , obj1_(obj1)
    , obj2_(obj2)
    , obj3_(obj3)
    , obj4_(obj4)
    , obj5_(obj5)
    , obj6_(obj6)
  { ENTER;
    obj1_.inc_ref();
    obj2_.inc_ref();
    obj3_.inc_ref();
    obj4_.inc_ref();
    obj5_.inc_ref();
    obj6_.inc_ref();
  }
  ~PyWrap6()
  { ENTER;
    obj1_.dec_ref();
    obj2_.dec_ref();
    obj3_.dec_ref();
    obj4_.dec_ref();
    obj5_.dec_ref();
    obj6_.dec_ref();
  }
 private:
  py::object obj1_;
  py::object obj2_;
  py::object obj3_;
  py::object obj4_;
  py::object obj5_;
  py::object obj6_;
};

template <typename T, typename W1, typename W2, typename W3, typename W4, typename W5, typename W6, typename W7>
class
__attribute__((visibility("hidden")))
PyWrap7 : public T
{
 public:
  PyWrap7(py::object obj1,
          py::object obj2,
          py::object obj3,
          py::object obj4,
          py::object obj5,
          py::object obj6,
          py::object obj7)
    : T(obj1.cast<W1*>(),
        obj2.cast<W2*>(),
        obj3.cast<W3*>(),
        obj4.cast<W4*>(),
        obj5.cast<W5*>(),
        obj6.cast<W6*>(),
        obj7.cast<W7*>())
    , obj1_(obj1)
    , obj2_(obj2)
    , obj3_(obj3)
    , obj4_(obj4)
    , obj5_(obj4)
    , obj6_(obj4)
    , obj7_(obj4)
  { ENTER;
    obj1_.inc_ref();
    obj2_.inc_ref();
    obj3_.inc_ref();
    obj4_.inc_ref();
    obj5_.inc_ref();
    obj6_.inc_ref();
    obj7_.inc_ref();
  }
  ~PyWrap7()
  { ENTER;
    obj1_.dec_ref();
    obj2_.dec_ref();
    obj3_.dec_ref();
    obj4_.dec_ref();
    obj5_.dec_ref();
    obj6_.dec_ref();
    obj7_.dec_ref();
  }
 private:
  py::object obj1_;
  py::object obj2_;
  py::object obj3_;
  py::object obj4_;
  py::object obj5_;
  py::object obj6_;
  py::object obj7_;
};

#endif /* PYWRAP_H_INCLUDED */
