// 2D Ray of polar coordinates and relationships with vector 2D operations,
// such as convertions, rotations, additions, etc.
//
#ifndef SYLVAIN__CODINGAME_RAY
#define SYLVAIN__CODINGAME_RAY

#include "vec2.hpp"

// The basic polar coordinate information, with an angle in degree and a radius
// equivalent to the distance to the pole. Radius is always positive. If it is
// found to be negative, program behaviour will be undefined.
struct Ray2 { int angle, rad; };

inline int angle(const Ray2& a) { return a.angle; }
inline int& angle(Ray2& a) { return a.angle; }

inline int rad(const Ray2& a) { return a.rad; }
inline int& rad(Ray2& a) { return a.rad; }

inline Ray2 operator+ (const Ray2& a, const Ray2& b)
{ return {angle(a) + angle(b), rad(a) + rad(b)}; }

inline Ray2 operator- (const Ray2& a, const Ray2& b)
{ return {angle(a) - angle(b), iabs(rad(a) - rad(b))}; }

inline Ray2 operator- (const Ray2& a)
{ return {-angle(a), rad(a)}; }

inline std::ostream& operator<< (std::ostream& o, const Ray2& a) {
  return o << "Ray2{" << angle(a) << ", " << rad(a) << "}";
}

inline Vec2 vec(const Ray2& a) { return Vec2{ icos(angle(a), rad(a)), isin(angle(a), rad(a)) }; }
inline Ray2 ray(const Vec2& a) { int r = mag(a); return Ray2{iacos3(x(a), y(a), r), r}; }

#endif // SYLVAIN__CODINGAME_RAY
