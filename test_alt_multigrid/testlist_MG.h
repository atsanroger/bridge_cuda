/*
        @file    $Id: testlist_MG.h #$

        @brief

        @author  Hideo Matsufuru
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-04-08 18:00:27 #$

        @version $LastChangedRevision: 2422 $
*/

#ifndef TESTLIST_MG_INCLUDED
#define TESTLIST_MG_INCLUDED

namespace test_sap_mult {
  int test(void);
}

//namespace test_afopr_coarse {
//  int test(void);
//  int test_parallel(void);
//}

namespace test_coarse {
  int test(void);
}

namespace test_restrict {
  int test(void);
}

namespace test_prolong {
  int test(void);
}

namespace test_MG_solver {
  int test(void);
}

#endif
