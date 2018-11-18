#define BOOST_TEST_MODULE Approximate Integer Tests
#include <boost/test/included/unit_test.hpp>

#include "math.hpp"
#include "ring.hpp"
#include "particle.hpp"

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

int error_isin_stdsin(int angle, int precision) {
  return iabs(isin(angle, precision) - int(std::sin(angle * M_PI / 180.0) * precision));
}

BOOST_AUTO_TEST_CASE(test_isin){
  constexpr const int SINE_PRECISION = 10000;
  BOOST_CHECK_EQUAL(isin(0, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(180, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(-180, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(360, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(-360, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(180, SINE_PRECISION), 0);
  BOOST_CHECK_EQUAL(isin(-180, SINE_PRECISION), 0);
  BOOST_CHECK_LT(error_isin_stdsin(30, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-30, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(60, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-60, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(90, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-90, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(120, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-120, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(150, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-150, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(210, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-210, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(240, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-240, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(270, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-270, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(300, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-300, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(330, SINE_PRECISION), 200);
  BOOST_CHECK_LT(error_isin_stdsin(-330, SINE_PRECISION), 200);
}

int iacos3_from_stdmath(int angle, int precision) {
  int x = std::cos(angle * M_PI / 180) * precision;
  int y = std::sin(angle * M_PI / 180) * precision;
  return iacos3(x, y, std::sqrt(x * x + y * y));
}

BOOST_AUTO_TEST_CASE(test_iacos2) {
  BOOST_CHECK_EQUAL(iacos2(1000, 0), 0);
  BOOST_CHECK_EQUAL(iacos2(-1000, 0), 180);
  BOOST_CHECK_EQUAL(iacos2(0, 1000), 90);
  BOOST_CHECK_EQUAL(iacos2(0, -1000), -90);
  BOOST_CHECK_LT(iabs(iacos2(500, 500)  - (  45)), 6);
  BOOST_CHECK_LT(iabs(iacos2(500, -500) - ( -45)), 6);
  BOOST_CHECK_LT(iabs(iacos2(-500, 500) - ( 135)), 6);
  BOOST_CHECK_LT(iabs(iacos3_from_stdmath(30, 10000))  - (  30), 10);
  BOOST_CHECK_LT(iabs(iacos3_from_stdmath(60, 10000))  - (  60), 10);
  BOOST_CHECK_LT(iabs(iacos3_from_stdmath(120, 10000)) - ( 120), 10);
  BOOST_CHECK_LT(iabs(iacos3_from_stdmath(150, 10000)) - ( 150), 10);
}

BOOST_AUTO_TEST_CASE(test_ring_anchor){
  Ring<int, 2> r(0);
  auto a = anchor<0>(r);
  auto b = anchor<1>(r);
  *a = 2;
  BOOST_CHECK_EQUAL(*a, 2);
  BOOST_CHECK_EQUAL(*b, 0);
  r.rotate();
  BOOST_CHECK_EQUAL(*a, 0);
  BOOST_CHECK_EQUAL(*b, 2);
}

template<typename Tp> void silence(const Tp&) { }

BOOST_AUTO_TEST_CASE(test_ring_move){
  Ring<int, 2> r1(2);
  BOOST_CHECK_EQUAL(r1.items()[0], 2);
  BOOST_CHECK_EQUAL(r1.items()[1], 2);
  {
    Ring<int, 2> r2(std::move(r1));
    silence(r2);
  }
  BOOST_CHECK_EQUAL(r1.items()[0], 0);
  BOOST_CHECK_EQUAL(r1.items()[1], 0);
}
