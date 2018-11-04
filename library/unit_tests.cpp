#define BOOST_TEST_MODULE Approximate Integer Tests
#include <boost/test/included/unit_test.hpp>

#include "math.hpp"

// BOOST pollutes my namespace! I'm thinking of raising a bug here.
BOOST_AUTO_TEST_CASE(test_iabs){
  BOOST_CHECK_EQUAL(iabs(0), 0);
  BOOST_CHECK_EQUAL(iabs(1), 1);
  BOOST_CHECK_EQUAL(iabs(-1), 1);
  BOOST_CHECK_EQUAL(iabs(1l), 1l);
  BOOST_CHECK_EQUAL(iabs(-1l), 1l);
}

BOOST_AUTO_TEST_CASE(test_irel){
  BOOST_CHECK_EQUAL(irel(0), 0);
  BOOST_CHECK_EQUAL(irel(1), 1);
  BOOST_CHECK_EQUAL(irel(-1), 0);
  BOOST_CHECK_EQUAL(irel(2l), 2l);
  BOOST_CHECK_EQUAL(irel(-2l), 0l);
}

BOOST_AUTO_TEST_CASE(test_amp){
  BOOST_CHECK_EQUAL(amp(0, 2), 2);
  BOOST_CHECK_EQUAL(amp(1, 2), 2);
  BOOST_CHECK_EQUAL(amp(-1, 2), 0);
  BOOST_CHECK_EQUAL(amp(2l, 2l), 2l);
  BOOST_CHECK_EQUAL(amp(-2l, 2l), 0l);
}

BOOST_AUTO_TEST_CASE(test_namp){
  BOOST_CHECK_EQUAL(namp(0, 2), 0);
  BOOST_CHECK_EQUAL(namp(1, 2), 0);
  BOOST_CHECK_EQUAL(namp(-1, 2), 2);
  BOOST_CHECK_EQUAL(namp(2l, 2l), 0l);
  BOOST_CHECK_EQUAL(namp(-2l, 2l), 2l);
}

BOOST_AUTO_TEST_CASE(test_ihyp){
  // Verifies that the error does not grow too much
  BOOST_CHECK_LT(ihyp(10, 10) - 14, 2);
  BOOST_CHECK_LT(ihyp(100, 100) - 141, 2);
  BOOST_CHECK_LT(ihyp(1000, 1000) - 1414, 4);
  BOOST_CHECK_LT(ihyp(10000, 10000) - 14142, 8);
  // Will overflow after
}

constexpr const int SINE_PRECISION = 10000;

BOOST_AUTO_TEST_CASE(test_isin){
  BOOST_CHECK_EQUAL(isin(0, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(90, 100), 1 * 100);
  BOOST_CHECK_EQUAL(isin(-90, SINE_PRECISION), -1 * SINE_PRECISION);
  BOOST_CHECK_EQUAL(isin(180, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(-180, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(270, SINE_PRECISION), -1 * SINE_PRECISION);
  BOOST_CHECK_EQUAL(isin(-270, SINE_PRECISION), 1 * SINE_PRECISION);
  BOOST_CHECK_EQUAL(isin(360, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(-360, SINE_PRECISION), 0);
  BOOST_CHECK_LT(iabs(isin(30, SINE_PRECISION) - int(std::sin(30.0 * M_PI / 180.0) * SINE_PRECISION)), 50);
  BOOST_CHECK_LT(iabs(isin(-30, SINE_PRECISION) - int(std::sin(-30.0 * M_PI / 180.0) * SINE_PRECISION)), 50);
  BOOST_CHECK_LT(iabs(isin(120, SINE_PRECISION) - int(std::sin(120.0 * M_PI / 180.0) * SINE_PRECISION)), 50);
  BOOST_CHECK_LT(iabs(isin(-120, SINE_PRECISION) - int(std::sin(-120.0 * M_PI / 180.0) * SINE_PRECISION)), 50);
  BOOST_CHECK_LT(iabs(isin(210, SINE_PRECISION) - int(std::sin(210.0 * M_PI / 180.0) * SINE_PRECISION)), 50);
  BOOST_CHECK_LT(iabs(isin(-210, SINE_PRECISION) - int(std::sin(-210.0 * M_PI / 180.0) * SINE_PRECISION)), 50);
  BOOST_CHECK_LT(iabs(isin(300, SINE_PRECISION) - int(std::sin(300.0 * M_PI / 180.0) * SINE_PRECISION)), 50);
  BOOST_CHECK_LT(iabs(isin(-300, SINE_PRECISION) - int(std::sin(-300.0 * M_PI / 180.0) * SINE_PRECISION)), 50);
}
