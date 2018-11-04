// 2-D Vector calculations and related functions.
//
// Specifically written for approximate integer calculations. I'm not sure why I
// made most of them as templates, cause `int` is what I use all the time. I
// guess I just like the opportunity cost of templates. :D
//

#include <array>
#include <iostream>
#include <cmath>

#include "math.hpp"

using namespace std;

struct Vec2 : array<int, 2>
{
    typedef array<int, 2> Base;

    Vec2() {}
    Vec2(int a, int b)
    { Base::operator[](0) = a; Base::operator[](1) = b; }

    int x() const { return Base::operator[](0); }
    int& x() { return Base::operator[](0); }

    int y() const { return Base::operator[](1); }
    int& y() { return Base::operator[](1); }
};

inline ostream& operator<< (ostream& o, const Vec2& a) {
  return o << "Vec2{" << a.x() << ", " << a.y() << "}";
}

inline Vec2 operator+ (const Vec2& a, const Vec2& b)
{ return {a.x() + b.x(), a.y() + b.y()}; }

inline Vec2 operator- (const Vec2& a, const Vec2& b)
{ return {a.x() - b.x(), a.y() - b.y()}; }

inline Vec2 operator- (const Vec2& a)
{ return {-a.x(), -a.y()}; }

inline Vec2 operator* (const Vec2& a, int factor)
{ return {a.x() * factor, a.y() * factor}; }

inline Vec2 operator/ (const Vec2& a, int factor)
{ return {a.x() / factor, a.y() / factor}; }

inline int distsq(const Vec2& a, const Vec2& b)
{ return sq(a.x() - b.x()) + sq(a.y() - b.y()); }

inline int magsq(const Vec2& a)
{ return sq(a.x()) + sq(a.y()); }

// Square root can't be computed on integers, but there's a fast convergent
// approximation with few iterations only: Newton's method.
//
//   - start at Manhattan magnitude (a.k.a taxicab metric)
//   - we know values are necessarily positive, so add +1 to avoid div by 0
//   - do 3 iterations of Newton's method
inline int mag3(const Vec2& a)
{
  int S = magsq(a);
  int x = abs(a.x()) + abs(a.y());
  x = (sq(x) + S) / (2 * x + 1);
  x = (sq(x) + S) / (2 * x + 1);
  return (sq(x) + S) / (2 * x + 1);
}

// Similar to mag3, normalize cares to:
//
//   - avoid overflows/underflows by computing divisions last
//   - avoid divisions by 0 with +1 since all values in the diviser are
//     guaranteed positives
inline Vec2 norm3(const Vec2& a, int norm)
{
  int x = mag3(a);
  return {(a.x() * norm) / (x + 1), (a.y() * norm) / (x + 1)};
}
