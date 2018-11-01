#include <array>

using namespace std;

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
