/*
        @file    testlist_accel.h
        @brief
        @author  Hideo Matsufuru
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2013-04-08 18:00:27 #$
        @version $LastChangedRevision: 2503 $
*/

#ifndef TESTLIST_ALT_INCLUDED
#define TESTLIST_ALT_INCLUDED


namespace Test_alt_Corelib {
  int test_all();
}

#ifdef USE_ALT_SIMD2
namespace Test_alt_SIMD2 {
  int test_all();
}
#endif

#ifdef USE_ALT_ACCEL
namespace Test_alt_Accel {
  int test_all();
}
#endif

#ifdef USE_ALT_VECTOR
namespace Test_alt_Vector {
  int test_all();
}
#endif

#ifdef USE_ALT_QXS
namespace Test_alt_QXS {
  int test_all();
}
#endif


#endif
