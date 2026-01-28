/*
        @file    $Id: testlist.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-04-08 18:00:27 #$

        @version $LastChangedRevision: 2528 $
*/

#ifndef TESTLIST_INCLUDED
#define TESTLIST_INCLUDED

//- prototype declarations


namespace Test_HMC_Staggered_Leapfrog {
  //  int update(void);
  int update_alt(void);
}

namespace Test_HMC_Domainwall {
  int update_Nf2_PV(void);
  int update_Nf2_PV_nosmr(void);
  int update_Nf2_PV_alt(void);
  int update_Nf2_PV_alt_Field(void);
  int update_Nf2_PV_alt_nosmr(void);
  int update_Nf2_PV_alt_eo_nosmr(void);
  int update_Nf2_PV_alt_eo(void);
  int update_Nf2p1_PV_alt(void);
  int update_Nf2p1_PV_alt_eo(void);

}

namespace Test_HMC_Wilson {
  int update_Nf2(void);
  int update_Nf2_nosmr(void);
  int update_Nf2_alt_Field(void);

  int update_Nf2_alt(void);
  int update_Nf2_alt_nosmr(void);

  int update_Nf2_alt_eo(void);
  int update_Nf2_alt_eo_ratio(void);
  int update_Nf2_alt_eo_nosmr(void);

  int RHMC_Nf2p1(void);
  int RHMC_Nf2p1_alt(void);
  int RHMC_Nf2p1_alt_eo(void);
  int RHMC_Nf2p1_alt_eo_ratio(void);
}

/*
namespace Test_HMC_Wilson {

}

namespace Test_HMC_Clover {
//  int update_Nf2(void);
  int run_test(void);
  int run_test_HYP(void);
  int update_Nf2_eo(void);
  int RHMC_Nf2p1(void);
  int RHMC_Nf2p1_2(void);  // HM added for test
  int RHMC_Nf2p1_3(void);  // HM added for test
  int RHMC_Nf2p1_4(void);  // HM added for test
  int RHMC_Nf2p1_5(void);  // HM added for test
  int RHMC_Nf2p1_eo(void);
  int RHMC_Nf2p1_eo_2(void);
  int RHMC_Nf2p1_eo_3(void);
}
*/

#endif
