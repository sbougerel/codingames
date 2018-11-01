// 2-D Vector calculations and related functions.
//
// Specifically written for approximate integer calculations. I'm not sure why I
// made most of them as templates, cause `int` is what I use all the time. I
// guess I just like the opportunity cost of templates. :D
//

#include <array>
#include <iostream>
#include <cmath>

using namespace std;

template<typename Tp>
constexpr Tp sq(Tp x) { return x * x; }

template<typename Tp>
struct Vec2 : array<Tp, 2>
{
    typedef array<Tp, 2> Base;

    Vec2() {}
    Vec2(Tp a, Tp b)
    { Base::operator[](0) = a; Base::operator[](1) = b; }

    Tp x() const { return Base::operator[](0); }
    Tp& x() { return Base::operator[](0); }

    Tp y() const { return Base::operator[](1); }
    Tp& y() { return Base::operator[](1); }
};

template<typename Tp>
inline ostream& operator<< (ostream& o, const Vec2<Tp>& a) {
  return o << "Vec2{" << a.x() << ", " << a.y() << "}";
}

template<typename Tp>
inline Vec2<Tp> operator+ (const Vec2<Tp>& a, const Vec2<Tp>& b)
{ return {a.x() + b.x(), a.y() + b.y()}; }

template<typename Tp>
inline Vec2<Tp> operator- (const Vec2<Tp>& a, const Vec2<Tp>& b)
{ return {a.x() - b.x(), a.y() - b.y()}; }

template<typename Tp>
inline Vec2<Tp> operator- (const Vec2<Tp>& a)
{ return {-a.x(), -a.y()}; }

template<typename Tp>
inline Vec2<Tp> operator* (const Vec2<Tp>& a, Tp factor)
{ return {a.x() * factor, a.y() * factor}; }

template<typename Tp>
inline Vec2<Tp> operator/ (const Vec2<Tp>& a, Tp factor)
{ return {a.x() / factor, a.y() / factor}; }

template<typename Tp>
inline Tp distsq(const Vec2<Tp>& a, const Vec2<Tp>& b)
{ return sq(a.x() - b.x()) + sq(a.y() - b.y()); }

template<typename Tp>
inline Tp magsq(const Vec2<Tp>& a)
{ return sq(a.x()) + sq(a.y()); }

// Square root can't be computed on integers, but there's a fast convergent
// approximation with few iterations only: Newton's method.
//
//   - start at Manhattan magnitude (a.k.a taxicab metric)
//   - we know values are necessarily positive, so add +1 to avoid div by 0
//   - do 3 iterations of Newton's method
template<typename Tp>
inline Tp mag3(const Vec2<Tp>& a)
{
  Tp S = magsq(a);
  Tp x = abs(a.x()) + abs(a.y());
  x = (sq(x) + S) / (2 * x + 1);
  x = (sq(x) + S) / (2 * x + 1);
  return (sq(x) + S) / (2 * x + 1);
}

// Similar to mag3, normalize cares to:
//
//   - avoid overflows/underflows by computing divisions last
//   - avoid divisions by 0 with +1 since all values in the diviser are
//     guaranteed positives
template<typename Tp>
inline Vec2<Tp> norm3(const Vec2<Tp>& a, Tp norm)
{
  Tp x = mag3(a);
  return {(a.x() * norm) / (x + 1), (a.y() * norm) / (x + 1)};
}
