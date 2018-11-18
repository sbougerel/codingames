#ifndef SYLVAIN__CODINGAME_PARTICLE
#define SYLVAIN__CODINGAME_PARTICLE

#include <iostream>

#include "vec2.hpp"
#include "ray2.hpp"

// The Particle object with the necessary trait accessors are defined below
//
struct Particle {
  Vec2 pos;
  Vec2 spd;
  int rad;
  int mass;
};

// The 4 pairs of accessor in the particule trait are defined below.  A particle
// as a position, a speed, a radius and a mass.
//
inline const Vec2& pos(const Particle& a) { return a.pos; }
inline Vec2& pos(Particle& a) { return a.pos; }
inline const Vec2& spd(const Particle& a) { return a.spd; }
inline Vec2& spd(Particle& a) { return a.spd; }
inline int rad(const Particle& a) { return a.rad; }
inline int& rad(Particle& a) { return a.rad; }
inline int mass(const Particle& a) { return a.mass; }
inline int& mass(Particle& a) { return a.mass; }

inline Particle free_move(const Particle& p) {
  Vec2 position = spd(p) + pos(p);
  return {position, spd(p), rad(p), mass(p)};
}

inline Particle free_move(const Particle& p, const Vec2& a) {
  Vec2 p_ = a / 2 + spd(p) + pos(p);
  Vec2 s_ = a + spd(p);
  return {p_, s_, rad(p), mass(p)};
}

// Policy for static dispatch during motion calculation. Needed since under
// constant policy we can do some nice optimizations in calcuation. And Constant
// acceleration is likely to happen very often.
struct ConstantThrustGradient { };
struct VariableThrustGradient { };

// Objects implementing the ThurstGradient trait only need to define at(), which
// must return a Vec2 object representing the acceleration to apply to a
// particule `p`. `iterations` is the number of steps to look into the future.
template<typename TG>
inline Vec2 at(const TG& tg, int iterations) {
  return tg.at(iterations);
}

template<typename TG>
inline Particle free_move(const Particle& p, const TG& t, int iterations, ConstantThrustGradient) {
  Vec2 a_ = at(t, 0);
  Vec2 p_ = (a_ * sq(iterations)) / 2 + spd(p) * iterations + pos(p);
  Vec2 s_ = a_ * iterations + spd(p);
  return {p_, s_, rad(p), mass(p)};
}

template<typename TG>
inline Particle free_move(const Particle& p, const TG& t, int iterations, VariableThrustGradient) {
  Vec2 p_ = pos(p);
  Vec2 s_ = spd(p);
  for (int i = 0; i < iterations; ++i) {
    Vec2 a_ = at(t, i);
    p_ = a_ / 2 + s_ + p_;
    s_ = a_ + s_;
  }
  return {p_, s_, rad(p), mass(p)};
}

template<typename TG>
inline Particle free_move(const Particle& p, const TG& t, int iterations) {
  // static dispatch
  return free_move(p, t, iterations, TG::ThrustGradiant());
}

// Zero thrust is here to be optimized away such that using this in a
// `free_move` function should compile to the same instruction as the
// `free_move` function without ThrustGradient;
struct ZeroThrust {
  typedef ConstantThrustGradient ThrustGradient;
  Vec2 at(int) const { return {0, 0}; }
};

// ConstantThrust is a ThurstGradiant implementation where the thrust will
// remain constant overtime.
struct ConstantThrust {
  typedef ConstantThrustGradient ThrustGradient;

  ConstantThrust(const Vec2& accel) : _accel(accel) { }
  Vec2 at(int) const { return _accel; }

private:
  Vec2 _accel;
};

// RotatingThrust is a ThurstGradiant implementation where the thrust will
// rotate overtime according to a constant angular speed (spin).
struct RotatingThrust {
  typedef VariableThrustGradient ThurstGratiant;

  RotatingThrust(const Vec2& thrust, int spin) : _thrust(ray(thrust)), _spin(spin) { }
  Vec2 at(int iterations) const
  { return vec(_thrust + Ray2{_spin * iterations, 0}); }

private:
  Ray2 _thrust;
  int _spin;
};

// Given two positions at discreet time t0 and t1 and a distance of closest
// appraoch `sqrad`, perform recursive halving of the time interval to check
// whether the particle collided.
//
inline int collide_int_sq(Vec2 x0, Vec2 y0, Vec2 x1, Vec2 y1, int sqrad) {
  constexpr const int stop_delta = 4;
  int sqd0 = distsq(x0, y0);
  if (sqd0 < sqrad) { return sqrad; }
  int sqd1 = distsq(x1, y1);
  while (true)
    {
      if (sqd1 < sqrad) { return sqrad; }
      if (   distsq(x0, x1) < stop_delta
          || distsq(y0, y1) < stop_delta) { break; }
      Vec2 xh = (x0 + x1) / 2;
      Vec2 yh = (y0 + y1) / 2;
      int sqd = distsq(xh, yh);
      if (sqd >= sqrad + magsq(x0 - xh) + magsq(y0 - yh)) { break; }
      if (sqd0 > sqd1)        // converge faster: chose the closest side
        { x0 = x1; y0 = y1; }
      x1 = xh; y1 = yh; sqd1 = sqd;
    }
  return (sqd0 < sqd1) ? sqd0 : sqd1;
}

// Collision Detection algorithm. Returns an a posteriori estimate of the
// closest approach between the 2 particles (squared) and the collision
// distance. When closest approach <= collision distance, there was
// collision.
//
// The estimation stops after max_iter, which is 100 by default, or when one of
// the particle gets out of the bounding box.
//
// The estimation is based on the ThrustGradiant associated with both
// particles.
//
template<typename TG1, typename TG2>
inline std::pair<int, int>
collide_sq (Particle x0, Particle y0, const TG1& tx, const TG2& ty,
            int max_iter = 100, const Box2& bb = {{-10000,-10000}, {10000, 10000}}) {
  int sqrad = sq(rad(x0)) + sq(rad(y0));
  int best_approach = distsq(pos(x0), pos(y0));
  if (best_approach <= sqrad) return std::make_pair(sqrad, sqrad);
  for (int i = 0; i < max_iter; ++ i) {
    Particle x1 = free_move(x0, at(tx, i));
    Particle y1 = free_move(y0, at(ty, i));
    // check against bb;
    int approach = collide_int_sq(pos(x0), pos(y0), pos(x1), pos(y1), sqrad);
    if (approach <= sqrad) return std::make_pair(sqrad, sqrad);
    if (approach < best_approach) { best_approach = approach; }
    x0 = x1; y0 = y1;
  }
  // return tuple with number of iteration at closest approach.
  return std::make_pair(sqrad, best_approach);
}

#endif // SYLVAIN__CODINGAME_PARTICLE
