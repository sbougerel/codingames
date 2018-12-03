#define BOOST_TEST_MODULE Approximate Integer Tests
#include <boost/test/included/unit_test.hpp>

#include "included.hpp"

// BOOST pollutes my namespace! I'm thinking of raising a bug here.
BOOST_AUTO_TEST_CASE(test_iabs){
  BOOST_CHECK_EQUAL(iabs(0), 0);
  BOOST_CHECK_EQUAL(iabs(1), 1);
  BOOST_CHECK_EQUAL(iabs(-1), 1);
  BOOST_CHECK_EQUAL(iabs(1), 1);
  BOOST_CHECK_EQUAL(iabs(-1), 1);
}

BOOST_AUTO_TEST_CASE(test_irel){
  BOOST_CHECK_EQUAL(irel(0), 0);
  BOOST_CHECK_EQUAL(irel(1), 1);
  BOOST_CHECK_EQUAL(irel(-1), 0);
  BOOST_CHECK_EQUAL(irel(2), 2);
  BOOST_CHECK_EQUAL(irel(-2), 0);
}

BOOST_AUTO_TEST_CASE(test_amp){
  BOOST_CHECK_EQUAL(amp(0, 2), 2);
  BOOST_CHECK_EQUAL(amp(1, 2), 2);
  BOOST_CHECK_EQUAL(amp(-1, 2), 0);
  BOOST_CHECK_EQUAL(amp(2, 2), 2);
  BOOST_CHECK_EQUAL(amp(-2, 2), 0);
}

BOOST_AUTO_TEST_CASE(test_namp){
  BOOST_CHECK_EQUAL(namp(0, 2), 0);
  BOOST_CHECK_EQUAL(namp(1, 2), 0);
  BOOST_CHECK_EQUAL(namp(-1, 2), 2);
  BOOST_CHECK_EQUAL(namp(2, 2), 0);
  BOOST_CHECK_EQUAL(namp(-2, 2), 2);
}

BOOST_AUTO_TEST_CASE(test_isgn){
  BOOST_CHECK_EQUAL(isgn(0, 0), 0);
  BOOST_CHECK_EQUAL(isgn(0, 2), 2);
  BOOST_CHECK_EQUAL(isgn(0, -2), -2);
  BOOST_CHECK_EQUAL(isgn(1, 0), 0);
  BOOST_CHECK_EQUAL(isgn(-1, 0), 0);
  BOOST_CHECK_EQUAL(isgn(1, 2), 2);
  BOOST_CHECK_EQUAL(isgn(-1, 2), -2);
  BOOST_CHECK_EQUAL(isgn(2, -2), -2);
  BOOST_CHECK_EQUAL(isgn(-2, -2), 2);
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
  BOOST_CHECK_LT(iabs(iacos2(500, 500)  - (  45)), 2);
  BOOST_CHECK_LT(iabs(iacos2(500, -500) - ( -45)), 2);
  BOOST_CHECK_LT(iabs(iacos2(-500, 500) - ( 135)), 2);
  BOOST_CHECK_LT(iabs(iacos3_from_stdmath(30, 10000))  - (  30), 2);
  BOOST_CHECK_LT(iabs(iacos3_from_stdmath(60, 10000))  - (  60), 2);
  BOOST_CHECK_LT(iabs(iacos3_from_stdmath(120, 10000)) - ( 120), 2);
  BOOST_CHECK_LT(iabs(iacos3_from_stdmath(150, 10000)) - ( 150), 2);
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

BOOST_AUTO_TEST_CASE(test_anorm){
  BOOST_CHECK_EQUAL(0,    anorm(0));
  BOOST_CHECK_EQUAL(0,    anorm(360));
  BOOST_CHECK_EQUAL(180,  anorm(180));
  BOOST_CHECK_EQUAL(-179, anorm(181));
  BOOST_CHECK_EQUAL(0,    anorm(-360));
  BOOST_CHECK_EQUAL(-180, anorm(-180));
  BOOST_CHECK_EQUAL(179,  anorm(-181));
  BOOST_CHECK_EQUAL(90,   anorm(90));
  BOOST_CHECK_EQUAL(-90,  anorm(-90));
  BOOST_CHECK_EQUAL(90,   anorm(450));
  BOOST_CHECK_EQUAL(-90,  anorm(-450));
}

BOOST_AUTO_TEST_CASE(test_adiff){
  BOOST_CHECK_EQUAL(0,    adiff(0, 0));
  BOOST_CHECK_EQUAL(0,    adiff(0, 360));
  BOOST_CHECK_EQUAL(0,    adiff(0, -360));
  BOOST_CHECK_EQUAL(0,    adiff(360, 0));
  BOOST_CHECK_EQUAL(0,    adiff(-360, 0));
  BOOST_CHECK_EQUAL(180,  adiff(180, 0));
  BOOST_CHECK_EQUAL(-180,  adiff(-180, 0));
  BOOST_CHECK_EQUAL(180,  adiff(0, -180));
  BOOST_CHECK_EQUAL(-180, adiff(0, 180));
}

BOOST_AUTO_TEST_CASE(test_adist){
  BOOST_CHECK_EQUAL(0,    adist(0, 0));
  BOOST_CHECK_EQUAL(0,    adist(0, 360));
  BOOST_CHECK_EQUAL(0,    adist(0, -360));
  BOOST_CHECK_EQUAL(0,    adist(360, 0));
  BOOST_CHECK_EQUAL(0,    adist(-360, 0));
  BOOST_CHECK_EQUAL(180,  adist(180, 0));
  BOOST_CHECK_EQUAL(180,  adist(-180, 0));
  BOOST_CHECK_EQUAL(180,  adist(0, -180));
  BOOST_CHECK_EQUAL(180,  adist(0, 180));
}

BOOST_AUTO_TEST_CASE(test_ray2_norm){
  BOOST_CHECK_EQUAL(Ray2({0, 0}),   norm(Ray2{0, 0}));
  BOOST_CHECK_EQUAL(Ray2({0, 0}),   norm(Ray2{360, 0}));
  BOOST_CHECK_EQUAL(Ray2({180, 0}),   norm(Ray2{180, 0}));
  BOOST_CHECK_EQUAL(Ray2({0, 0}),   norm(Ray2{-360, 0}));
  BOOST_CHECK_EQUAL(Ray2({-180, 0}),   norm(Ray2{-180, 0}));
  BOOST_CHECK_EQUAL(Ray2({90, 0}),  norm(Ray2{90, 0}));
  BOOST_CHECK_EQUAL(Ray2({-90, 0}), norm(Ray2{-90, 0}));
  BOOST_CHECK_EQUAL(Ray2({90, 0}),  norm(Ray2{450, 0}));
  BOOST_CHECK_EQUAL(Ray2({-90, 0}), norm(Ray2{-450, 0}));
}

BOOST_AUTO_TEST_CASE(test_free_move) {
  // Make 2 particles at the edge of the board, facing each other
  Particle x0 = {{-10000, 0}, {100, 0}, 0, 500, 1};
  Physics<RealisticThrustModel, VaccumDragModel> model;
  // The second implementation is actually less precise.
  BOOST_CHECK_LT(iabs(x(pos(reaction(x0, {0, 0}, 100)))
                      - x(pos(iterate_reaction(100, x0, CoastingAction(), model)))),
                 100);
}

BOOST_AUTO_TEST_CASE(test_collide_int) {
  // At bounds
  BOOST_CHECK_EQUAL(linear_collide(Vec2({0, 0}), Vec2({0, 0}),
                                   Vec2({0, 0}), Vec2({0, 0}), 0),
                    0);
  // Parallel, never really meets: check against forever loops
  BOOST_CHECK_EQUAL(linear_collide(Vec2({0, 0}), Vec2({0, 1000}),
                                   Vec2({1000, 0}), Vec2({1000, 1000}), 1000000),
                    1000000);
  // Face each others: collide
  BOOST_CHECK_EQUAL(linear_collide(Vec2({0, 0}), Vec2({0, 1000}),
                                   Vec2({0, 1000}), Vec2({0, 0}), 100),
                    100);
  // Cross each others: collide
  BOOST_CHECK_EQUAL(linear_collide(Vec2({-10000, 0}), Vec2({0, 10000}),
                                   Vec2({10000, 0}), Vec2({0, -10000}),
                                   100),
                    100);
  // Follow each other: miss
  BOOST_CHECK_EQUAL(linear_collide(Vec2({-10000, 0}), Vec2({0, 0}),
                                   Vec2({0, 0}), Vec2({10000, 0}),
                                   10000000),
                    100000000);
}

BOOST_AUTO_TEST_CASE(test_collide_two) {
  // Make 2 particles at the edge of the board, facing each other
  Particle x0 = {{-10000, 0}, {100, 0}, 0, 500, 1};
  Particle x1 = {{10000, 0}, {-100, 0}, 0, 500, 1};
  Physics<InstantThrustModel, VaccumDragModel> model1;
  Physics<RealisticThrustModel, VaccumDragModel> model2;
  // In vaccum, without acceleration, they should collide
  BOOST_CHECK_EQUAL(std::get<1>(collide_two(x0, x1, CoastingAction(), CoastingAction(), model1)),
                    sq(500) + sq(500));
  // Now with a constant acceleration, opposite to the speed, they get close but
  // do not touch. Compare with the similarly imprecise version.
  BOOST_CHECK_EQUAL(std::get<1>(collide_two(x0, x1, ConstantAction({-10, 0}), ConstantAction({10, 0}), model2, 10)),
                    distsq(pos(reaction(x0, {-10, 0}, 10)),
                           pos(reaction(x1, {10, 0}, 10))));
  // Now aiming for a collision at a lower point (still moving into each others)
  BOOST_CHECK_EQUAL(std::get<1>(collide_two(x0, x1, TargetAction({0, 2000}, 10), TargetAction({0, 2000}, 10), model1)),
                    sq(500) + sq(500));
  // Now aiming for avoidance
  BOOST_CHECK_NE(std::get<1>(collide_two(x0, x1, TargetAction({0, -2000}, 10), TargetAction({0, 2000}, 10), model1)),
                 sq(500) + sq(500));
  // Now reaching the final position due to exessive drag
  Physics<RealisticThrustModel, BasicDragModel<20, 300>> model3;
  BOOST_CHECK_NE(std::get<1>(collide_two(x0, x1, TargetAction({0, -2000}, 20), TargetAction({0, 2000}, 20), model3)),
                 sq(500) + sq(500));
}
