#ifndef SYLVAIN__CODINGAME_PARTICLE
#define SYLVAIN__CODINGAME_PARTICLE

#include <iostream>

#include "vec2.hpp"

// The 4 pairs of accessor in the particule trait are defined below.  A particle
// as a position, a speed, a radius and a mass. These 3 properties are
// sufficient to define - even when assuming arbitrary units - collisions,
// approach, etc.
//
template<typename Tp>
inline const Vec2& pos(const Tp& a) { return a.pos(); }
template<typename Tp>
inline Vec2& pos(Tp& a) { return a.pos(); }

template<typename Tp>
inline const Vec2& spd(const Tp& a) { return a.spd(); }
template<typename Tp>
inline Vec2& spd(Tp& a) { return a.spd(); }

template<typename Tp>
inline int rad(const Tp& a) { return a.rad(); }
template<typename Tp>
inline int& rad(Tp& a) { return a.rad(); }

template<typename Tp>
inline int mass(const Tp& a) { return a.rad(); }
template<typename Tp>
inline int& mass(Tp& a) { return a.rad(); }

template<typename Tp>
inline int spin(const Tp& a) { return a.spin(); }
template<typename Tp>
inline int& spin(Tp& a) { return a.spin(); }

// The Particle object with the necessary trait overload are defined below
//
struct Particle {
  explicit Particle(const Vec2& pos_ = Vec2(), const Vec2& spd_ = Vec2(),
                    int rad_ = 0, int spin_ = 0, int mass_ = 0)
    : pos(pos_), spd(spd_), rad(rad_), spin(spin_), mass(mass_) { }
  Vec2 pos;
  Vec2 spd;
  int rad;
  int spin;
  int mass;
};

inline const Vec2& pos(const Particle& a) { return a.pos; }
inline Vec2& pos(Particle& a) { return a.pos; }

inline const Vec2& spd(const Particle& a) { return a.spd; }
inline Vec2& spd(Particle& a) { return a.spd; }

inline int rad(const Particle& a) { return a.rad; }
inline int& rad(Particle& a) { return a.rad; }

inline int mass(const Particle& a) { return a.mass; }
inline int& mass(Particle& a) { return a.mass; }

inline int spin(const Particle& a) { return a.spin; }
inline int& spin(Particle& a) { return a.spin; }

// Linear Collision Detection algorithm. A posteriori (not exact). Returns 0 if
// the 2 particles collide, otherwise gives an a posteriori estimate of the
// closest approach between the 2 particles.
//
// The estimation is based on a constant speed for both particles x and y.
template<typename P>
inline int collide_straight(P x, P y) {
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
// the 2 particles collide, otherwise gives an a posteriori estimate of the
// closest approach between the 2 particles.
//
// The estimation is based on a constant speed *and a constant spin* for both
// particles x and y.
template<typename P>
inline int collide_curve(P x, P y) {
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
