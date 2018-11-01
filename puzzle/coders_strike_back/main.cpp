#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
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

template<typename Tp>
inline Tp mag3(const Vec2<Tp>& a)
{
  Tp S = magsq(a);
  Tp x = abs(a.x()) + abs(a.y());
  x = (sq(x) + S) / (2 * x + 1);
  x = (sq(x) + S) / (2 * x + 1);
  return (sq(x) + S) / (2 * x + 1);
}

template<typename Tp>
inline Vec2<Tp> norm3(const Vec2<Tp>& a, Tp norm)
{
  Tp x = mag3(a);
  return {(a.x() * norm) / (x + 1), (a.y() * norm) / (x + 1)};
}


/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/
int main()
{
    bool boost_used = false;

    // game loop
    while (1) {
        int x;
        int y;
        int nextCheckpointX; // x position of the next check point
        int nextCheckpointY; // y position of the next check point
        int nextCheckpointDist; // distance to the next checkpoint
        int nextCheckpointAngle; // angle between your pod orientation and the direction of the next checkpoint
        cin >> x >> y >> nextCheckpointX >> nextCheckpointY >> nextCheckpointDist >> nextCheckpointAngle; cin.ignore();
        int opponentX;
        int opponentY;
        cin >> opponentX >> opponentY; cin.ignore();

        int thrust;
        if (nextCheckpointAngle > 100 || nextCheckpointAngle < -100)
            cout << nextCheckpointX << " " << nextCheckpointY << " 0" << endl;
        else
            if (!boost_used && (nextCheckpointAngle < 10 && nextCheckpointAngle > -10) && nextCheckpointDist > 5000) {
                cout << nextCheckpointX << " " << nextCheckpointY << " BOOST" << endl;
                boost_used = true;
            }
            else
                cout << nextCheckpointX << " " << nextCheckpointY << " 100" << endl;
    }
}
