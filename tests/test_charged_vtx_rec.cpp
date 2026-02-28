#define BOOST_TEST_MODULE ChargedVtxRecTests
#include <boost/test/unit_test.hpp>
#include <cmath>

#include "../Include/Codes/inc/charged_mom.h"
#include "../Include/Codes/inc/const.h"

BOOST_AUTO_TEST_CASE(charged_mom_mode1_basic)
{
  ErrorHandling::ErrorLogs logger("/tmp/kloe_test_logs/", 0, "");
  KLOE::ChargedVtxRec<double, int> rec(logger);

  double mom[4] = {0.0, 0.0, 0.0, 0.0};
  rec.charged_mom(2.0, 0.0, 1.0, mom, 1, logger);

  BOOST_CHECK_CLOSE(mom[0], 500.0, 1e-6);
  BOOST_CHECK_SMALL(mom[1], 1e-9);
  BOOST_CHECK_CLOSE(mom[2], 500.0, 1e-6);

  const double p2 = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
  const double expectedE = std::sqrt(p2 + PhysicsConstants::mPiCh * PhysicsConstants::mPiCh);
  BOOST_CHECK_CLOSE(mom[3], expectedE, 1e-6);
}

BOOST_AUTO_TEST_CASE(ipboostcorr_denom_zero)
{
  ErrorHandling::ErrorLogs logger("/tmp/kloe_test_logs/", 0, "");
  KLOE::ChargedVtxRec<double, int> rec(logger);

  double X_line[3] = {0.0, 0.0, 0.0};
  double vec_line[3] = {1.0, 0.0, 0.0};
  double X_plane[3] = {0.0, 0.0, 0.0};
  double vec_plane[3] = {0.0, 1.0, 0.0};
  double out[3] = {0.0, 0.0, 0.0};

  const int rc = rec.IPBoostCorr(X_line, vec_line, X_plane, vec_plane, out);
  BOOST_CHECK_EQUAL(rc, static_cast<int>(ErrorHandling::ErrorCodes::DENOM_EQ_ZERO));
}