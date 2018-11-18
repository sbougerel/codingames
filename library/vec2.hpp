// 2-D Vector calculations and related functions.
//
// Specifically written for approximate integer calculations. I'm not sure why I
// made most of them as templates, cause `int` is what I use all the time. I
// guess I just like the opportunity cost of templates. :D
//
#ifndef SYLVAIN__CODINGAME_VEC2
#define SYLVAIN__CODINGAME_VEC2

#include <array>
#include <iostream>
#include <cmath>

#include "math.hpp"

struct Vec2 { int x, y; };

inline int x(const Vec2& a) { return a.x; }
inline int& x(Vec2& a) { return a.x; }

inline int y(const Vec2& a) { return a.y; }
inline int& y(Vec2& a) { return a.y; }

inline std::ostream& operator<< (std::ostream& o, const Vec2& a) {
  return o << "Vec2{" << x(a) << ", " << y(a) << "}";
}

inline Vec2 operator+ (const Vec2& a, const Vec2& b)
{ return {x(a) + x(b), y(a) + y(b)}; }

inline Vec2 operator- (const Vec2& a, const Vec2& b)
{ return {x(a) - x(b), y(a) - y(b)}; }

inline Vec2 operator- (const Vec2& a)
{ return {-x(a), -y(a)}; }

inline Vec2 operator* (const Vec2& a, int factor)
{ return {x(a) * factor, y(a) * factor}; }

inline Vec2 operator/ (const Vec2& a, int factor)
{ return {x(a) / factor, y(a) / factor}; }

inline Vec2 operator<< (const Vec2& a, int factor)
{ return {x(a) << factor, y(a) << factor}; }

inline Vec2 operator>> (const Vec2& a, int factor)
{ return {x(a) >> factor, y(a) >> factor}; }

inline int magsq(const Vec2& a)
{ return sq(x(a)) + sq(y(a)); }

inline int distsq(const Vec2& a, const Vec2& b)
{ return magsq(a - b); }

// Square root can't be computed on integers, but there's a fast convergent
// approximation with few iterations only: Newton's method.
//
//   - start at Manhattan magnitude (a.k.a taxicab metric)
//   - we know values are necessarily positive, so add +1 to avoid div by 0
//   - do 3 iterations of Newton's method
inline int mag(const Vec2& a)
{
  int S = magsq(a);
  int m = iabs(x(a)) + iabs(y(a));
  m = (sq(m) + S) / (2 * m + 1);
  m = (sq(m) + S) / (2 * m + 1);
  return (sq(m) + S) / (2 * m + 1);
}

// Similar to mag, normalize cares to:
//
//   - avoid overflows/underflows by computing divisions last
//   - avoid divisions by 0 with +1 since all values in the diviser are
//     guaranteed positives
inline Vec2 norm(const Vec2& a, int norm)
{
  int m = mag(a);
  return {(x(a) * norm) / (m + 1), (y(a) * norm) / (m + 1)};
}

#endif // SYLVAIN__CODINGAME_VEC2
