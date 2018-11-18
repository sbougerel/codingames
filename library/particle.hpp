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

inline Particle free_move(const Particle& p, int iteration) {
  Vec2 position = spd(p) * iteration + pos(p);
  return {position, spd(p), rad(p), mass(p)};
}

// Policy for static dispatch during motion calculation. Needed since under
// constant policy we can do some nice optimizations in calcuation. And Constant
// acceleration is likely to happen very often.
struct ConstantThrustGradient { };
struct VariableThrustGradient { };

// Objects implementing the ThurstGradient trait only need to define at(), which
// must return a Vec2 object representing the acceleration to apply to a
// particule `p`. `iterations` is the number of steps to look into the future.
template<typename TGImpl>
inline Vec2 at(const TGImpl& tg, int iterations) {
  return tg.at(iterations);
}

template<typename TG>
inline Particle free_move(const Particle& p, const TG& t) {
  Vec2 a_ = at(t, 0);
  Vec2 p_ = a_ / 2 + spd(p) + pos(p);
  Vec2 s_ = a_ + spd(p);
  return {p_, s_, rad(p), mass(p)};
}

template<typename TG>
inline Particle free_move(const Particle& p, const TG& t, int iterations, ConstantThrustGradient) {
  Vec2 a_ = at(t, 0);
  Vec2 p_ = (a_ * sq(iterations)) / 2 + spd(p) * iterations + pos(p);
  Vec2 s_ = a_ * iterations + spd(p);
  return {p_, s_, rad(p), mass(p)};
}

template<typename P, typename T>
inline P free_move(const P& p, const T& t, int iterations, VariableThrustGradient) {
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
// rotate over time according to a constant angular speed (spin).
struct RotatingThrust {
  typedef VariableThrustGradient ThurstGratiant;

  RotatingThrust(const Vec2& thrust, int spin) : _thrust(ray(thrust)), _spin(spin) { }
  Vec2 at(int iterations) const
  { return vec(_thrust + Ray2{_spin * iterations, 0}); }

private:
  Ray2 _thrust;
  int _spin;
};

// Linear Collision Detection algorithm. A posteriori (not exact). Returns 0 if
// the 2 particles collide, otherwise gives an a posteriori estimate of the
// closest approach between the 2 particles.
//
// The estimation is based on a constant speed for both particles x and y.
inline int collide_straight(Particle x, Particle y) {
  // iterate so long as the particules are coming closer to each others
  int prev_dist;
  int curr_dist = distsq(pos(x), pos(y));
  do
    {
      prev_dist = curr_dist;
      Vec2 prev_x = pos(x);
      Vec2 prev_y = pos(y);
      pos(x) = pos(x) + spd(x);
      pos(y) = pos(y) + spd(y);
      curr_dist = distsq(pos(x), pos(y));
      if (curr_dist < sq(rad(x)) + sq(rad(y))) { return 0; }
      else {
        // check collision between virtual particles spanning the entire motion
        Vec2 spd_x = spd(x) >> 1;
        Vec2 spd_y = spd(y) >> 1;
        Vec2 pos_x = prev_x + spd_x;
        Vec2 pos_y = prev_y + spd_y;
        int dist = distsq(pos_x, pos_y);
        if (dist < sq(rad(x)) + sq(rad(y)) + magsq(spd_x) + magsq(spd_y)) {
          spd(x) = spd_x;
          spd(y) = spd_y;
          pos(x) = pos_x;
          pos(y) = pos_y;
          curr_dist = dist;
        }
      }
    } while (prev_dist > curr_dist);
  return prev_dist;
}

// Linear Collision Detection algorithm. A posteriori (not exact). Returns 0 if
// the two particles collide, otherwise gives an a posteriori estimate of the
// closest approach between the 2 particles.
//
// The estimation is based on a constant acceleration *and a constant spin* for both
// particles x and y.
template<typename P>
inline int collide(P x, P y, int max_iter = 1000) {
  // iterate so long as the particules are coming closer to each others
  int prev_dist;
  int curr_dist = distsq(pos(x), pos(y));
  do
    {
      prev_dist = curr_dist;
      Vec2 prev_x = pos(x);
      Vec2 prev_y = pos(y);
      pos(x) = pos(x) + spd(x);
      pos(y) = pos(y) + spd(y);
      curr_dist = distsq(pos(x), pos(y));
      if (curr_dist < sq(rad(x)) + sq(rad(y))) { return 0; }
      else {
        // check collision between virtual particles spanning the entire motion
        Vec2 spd_x = spd(x) >> 1;
        Vec2 spd_y = spd(y) >> 1;
        Vec2 pos_x = prev_x + spd_x;
        Vec2 pos_y = prev_y + spd_y;
        int dist = distsq(pos_x, pos_y);
        if (dist < sq(rad(x)) + sq(rad(y)) + magsq(spd_x) + magsq(spd_y)) {
          spd(x) = spd_x;
          spd(y) = spd_y;
          pos(x) = pos_x;
          pos(y) = pos_y;
          curr_dist = dist;
        }
      }
    } while (prev_dist > curr_dist);
  return prev_dist;
}

#endif // SYLVAIN__CODINGAME_PARTICLE
