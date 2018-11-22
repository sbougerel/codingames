#ifndef SYLVAIN__CODINGAME_INCLUDED
#define SYLVAIN__CODINGAME_INCLUDED

#include <tuple>
#include <cmath>
#include <array>
#include <iostream>

// This is the canonical implementation for absolute. It's so beautiful, I
// wanted to write it again. It just blows my mind everytime I look at it.
//
// But the std library abs() is even smarter:
//
//      <+0>:   push   %rbp
//      <+1>:   mov    %rsp,%rbp
//      <+4>:   mov    %rdi,-0x8(%rbp)
//      <+8>:   mov    -0x8(%rbp),%rax
//     <+12>:   cqto
//     <+14>:   mov    %rdx,%rax
//     <+17>:   xor    -0x8(%rbp),%rax
//     <+21>:   sub    %rdx,%rax
//     <+24>:   pop    %rbp
//     <+25>:   retq
//
// It uses `cqto` to convert long to double long, and boom: it also gets the
// mask in another register in a single operation 👏!
inline int iabs(int a) {
  int t = a >> (sizeof(int) * 8 - 1);
  return (a^t) - t;
}

// irel() acts as linear rectifier; it returns identity for positive integers,
// but returns 0 for negative integers, such as:
//
//     if (x > 0) return x;
//     else return 0;
inline int irel(int a) {
  int t = a >> (sizeof(int) * 8 - 1);
  return (a|t) - t;
}

// nirel() returns identity for negative integers, but returns 0 for positive
// integers:
//
//     if (x < 0) return x;
//     else return 0;
inline int nirel(int a) {
  int t = a >> (sizeof(int) * 8 - 1);
  return (a & t);
}

// isgn() returns `boost` if `gate` is positive or zero, and `-boost` if `gate`
// is negative for integers:
//
//     if (gate >= 0) return boost;
//     if (gate < 0) return -boost;
inline int isgn(int gate, int boost) {
  int t = gate >> (sizeof(int) * 8 - 1);
  return (boost|t) - t + (t & (-boost));
}

// isgv() returns `pos` if `gate` is positive or zero, and `neg` if `gate` is
// strictly negative for integers:
//
//     if (gate >= 0) return pos;
//     if (gate < 0) return neg;
inline int isgv(int gate, int pos, int neg) {
  int t = gate >> (sizeof(int) * 8 - 1);
  return (pos|t) - t + (t & neg);
}

// amp() returns `boost` if `gate` is a positive integer or 0, and returns
// `boost` otherwise for integers in complement 2 notation;
//
//     if (gate >= 0) return boost;
//     else return 0;
inline int amp(int gate, int boost) {
  int t = gate >> (sizeof(int) * 8 - 1);
  return (boost|t) - t;
}

// namp() returns `boost` if `gate` is a striclty negative integer, and 0
// otherwise for integers in complement 2 notation;
//
//     if (x < 0) return boost;
//     else return 0;
inline int namp(int gate, int boost) {
  int t = gate >> (sizeof(int) * 8 - 1);
  return (t & boost);
}

// Squares a value.
template<typename Tp>
constexpr Tp sq(Tp x) { return x * x; }

// Approximate hypotenuse implementation. Square root can't be computed on
// integers, but there's a fast convergent approximation with few iterations
// only: Newton's method.
//
//   - start at the Manhattan magnitude (a.k.a taxicab metric)
//   - we know values are necessarily positive, so add +1 to avoid div by 0
//   - do 3 iterations of Newton's method
inline int ihyp(const int adjacent, const int opposite)
{
  int S = sq(adjacent) + sq(opposite);
  int x = iabs(adjacent) + iabs(opposite);
  x = (sq(x) + S) / (2 * x + 1);
  x = (sq(x) + S) / (2 * x + 1);
  return (sq(x) + S) / (2 * x + 1);
}

// Return the sine value for an `angle` in degree, multipled by arbitrary
// precision.
//
// It does it by mapping degree angles into into 512 values. The first 7 bits
// are the angle within a 90 degree quadrant. The next 2 bits control the choice
// of sign and quadrant to use: sin(x) or sin(90 - x).
//
// It should have less than 2% of error.
inline int isin(int angle, int scale) {
  constexpr const int Factor = 81;               // 256 / PI ~= 81
  constexpr const int Factor2 = 54 * sq(Factor); // 9 * 3! * (256/PI ~= 81)^2
  constexpr const float Factor3 = Factor2 * Factor;
  constexpr const int bits = sizeof(int) * 8;
  int s = angle >> (bits - 1);                   // sign mask: 0 or -1
  int aa = (((angle^s) - s) * 128 + 45) / 90;    // absolute angle in 256-unit per PI
  int h = (aa << (bits - 9)) >> (bits - 1);      // first or second half: 0 or -1
  int ra = aa & 0x7F;                            // 128-unit angle in quadrant
  if (aa & 0x80) { ra = 128 - ra; }              // first or second quadrant
  float f = float(ra * (Factor2 - (sq(ra) * 8))) / Factor3;
  ra = scale * f;
  s = h^s;                                       // XOR halves and sign together
  return (ra^s) - s;
}

inline int icos(int angle, int scale) {
  return isin(90 - angle, scale);
}

// Return an angle in degree for the value of the adjacent length `x`, the
// opposite length`y` and their hypotenuse (when already precomputed), with less
// than 3% error (10 degrees over 360 degrees).
//
// Uses the lagrange approximation that can be found everywhere on Internet,
// only adapted to integers and degrees. The input range is only defined between
// [0, 1], which means we need to use the right terms for x vs y:
//
//     (-0.698131700 * x * x -0.872664625) * x + 1.570796325;
//
// if hypot is 0, the result is undefined.
inline int iacos3(int x, int y, int hypot) {
  constexpr const float A = -0.698131700f;
  constexpr const float B = -0.872664625f;
  constexpr const float C = A + B;
  constexpr const float D = 180.f / M_PI;
  float f = float(x) / float(hypot);
  f = (A * f * f + B) * f - C;
  int r = (D * f);
  return isgn(y, r);
}

// Return an angle in degree for the value of the adjacent length `x`, the
// opposite length`y` with less than 3% error.
//
inline int iacos2(int x, int y) {
  return iacos3(x, y, ihyp(x, y));
}

struct Vec2 { int x, y; };

inline int x(const Vec2& a) { return a.x; }
inline int& x(Vec2& a) { return a.x; }

inline int y(const Vec2& a) { return a.y; }
inline int& y(Vec2& a) { return a.y; }

inline std::ostream& operator<< (std::ostream& o, const Vec2& a) {
  return o << "Vec2({" << x(a) << ", " << y(a) << "})";
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

inline bool operator== (const Vec2& a, const Vec2& b)
{ return (x(a) == x(b) && y(a) == y(b)); }

inline bool operator!= (const Vec2& a, const Vec2& b) { return !(a == b); }

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

struct Box2 { Vec2 low, high; };

inline const Vec2& low(const Box2& b) { return b.low; }
inline const Vec2& high(const Box2& b) { return b.high; }
inline Vec2& low(Box2& b) { return b.low; }
inline Vec2& high(Box2& b) { return b.high; }

inline bool within(const Box2& b, const Vec2& v) {
  return (x(v) >= x(low(b)) && x(v) < x(high(b))
          && y(v) >= y(low(b)) && y(v) < y(high(b)));
}

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

inline bool operator== (const Ray2& a, const Ray2& b)
{ return (angle(a) == angle(b) && rad(a) == rad(b)); }

inline bool operator!= (const Ray2& a, const Ray2& b) { return !(a == b); }

inline std::ostream& operator<< (std::ostream& o, const Ray2& a) {
  return o << "Ray2({" << angle(a) << ", " << rad(a) << "})";
}

inline Ray2 norm(const Ray2& a) { return {angle(a) % 360, rad(a)}; }

inline Vec2 vec(const Ray2& a) { return Vec2{ icos(angle(a), rad(a)), isin(angle(a), rad(a)) }; }
inline Ray2 ray(const Vec2& a) {
  int r = mag(a);
  return (r == 0) ? Ray2{0, 0} : Ray2{iacos3(x(a), y(a), r), r};
}

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

inline bool operator==(const Particle& a, const Particle& b) {
  return (pos(a) == pos(b)
          && spd(a) == spd(b)
          && rad(a) == rad(b)
          && mass(a) == mass(b));
}
inline bool operator!=(const Particle& a, const Particle& b)
{ return !(a == b); }

inline std::ostream& operator<<(std::ostream& o, const Particle& a) {
  o << "Particle({" << pos(a) << ", " << spd(a) << ", " << rad(a) << ", " << mass(a) << "})";
  return o;
}

inline Particle free_move(const Particle& p) {
  return {spd(p) + pos(p), spd(p), rad(p), mass(p)};
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
  return free_move(p, t, iterations, typename TG::ThrustGradient());
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
  typedef VariableThrustGradient ThrustGradient;

  RotatingThrust(const Ray2& thrust, int spin) : _thrust(thrust), _spin(spin) { }
  Vec2 at(int iterations) const
  { return vec(Ray2{_spin * iterations + angle(_thrust), rad(_thrust)}); }

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

// Collision Detection algorithm. Returns the distance of collision, an a
// posteriori estimate of the closest approach between the 2 particles (squared)
// and the time to collision. When closest approach <= distance of collision,
// collision has occured.
//
// The estimation stops after max_iter, which is 100 by default, or when one of
// the particle gets out of the bounding box.
//
// The estimation is based on the ThrustGradiant associated with both
// particles.
//
template<typename TG1, typename TG2>
inline std::tuple<int, int, int>
collide_sq(Particle x0, Particle y0, const TG1& tx, const TG2& ty,
            int max_iter = 100, const Box2& bb = {{-10000,-10000}, {10000, 10000}}) {
  int sqrad = sq(rad(x0)) + sq(rad(y0));
  int best_approach = distsq(pos(x0), pos(y0));
  if (best_approach <= sqrad)
    { return std::make_tuple(sqrad, best_approach, 0); }
  int i = 0;
  for (; i < max_iter; ++i) {
    Particle x1 = free_move(x0, at(tx, i));
    Particle y1 = free_move(y0, at(ty, i));
    if (!within(bb, pos(x1)) || !within(bb, pos(y1))) break;
    int approach = collide_int_sq(pos(x0), pos(y0), pos(x1), pos(y1), sqrad);
    if (approach <= sqrad)
      { return std::make_tuple(sqrad, approach, i); }
    if (approach < best_approach) { best_approach = approach; }
    x0 = x1; y0 = y1;
  }
  return std::make_tuple(sqrad, best_approach, i);
}

// Ring & Anchor are 2 simple objects that store objects in a contiguous location
// and then rotate addressing to each objects stored.
//
// Ring is designed to store a small amount of large types. It is non-copyable
// but movable.
template<typename Tp, unsigned N>
struct Ring {
  Ring() { _initialize_aquire(); }
  explicit Ring(const Tp& value) { _initialize_aquire(value); }
  Ring(const Ring& ring) = delete;
  Ring(Ring&& ring) {
    _initialize_aquire();
    std::swap(_data, ring._data);
    for (unsigned i = 0; i < N; ++i) { std::swap(_addr[i], ring._addr[i]); }
  }
  ~Ring() { delete _data; }

  void rotate() {
    for (unsigned i = N - 1; i > 0; --i) { std::swap(_addr[i], _addr[i - 1]); }
  }

  std::array<Tp, N>& items() { return *_data; }
  const std::array<Tp, N>& items() const { return *_data; }

  std::array<Tp*, N>& addresses() { return _addr; }
  const std::array<Tp*, N>& addresses() const { return _addr; }

  template<unsigned P>
  struct Anchor {
    Anchor(Ring<Tp, N>& r) : _ref(r) { }
    Tp& operator*() { return *std::get<P>(_ref.addresses()); }
    Tp* operator->() { return std::get<P>(_ref.addresses()); }

  private:
    Ring<Tp, N>& _ref;
  };

  template<unsigned P>
  struct ConstAnchor {
    ConstAnchor(const Ring<Tp, N>& r) : _ref(r) { }
    const Tp& operator*() { return *std::get<P>(_ref.addresses()); }
    const Tp* operator->() { return std::get<P>(_ref.addresses()); }

  private:
    const Ring<Tp, N>& _ref;
  };

private:
  void _initialize_aquire(const Tp& val = Tp()) {
    _data = new std::array<Tp, N>();
    for (unsigned i = 0; i < N; ++i) { _addr[i] = &(*_data)[i]; (*_data)[i] = val; }
  }
  std::array<Tp*, N> _addr;
  std::array<Tp, N>* _data;
};

// Class undefined for 0 elements
template<typename Tp> struct Ring<Tp, 0> { };

template<unsigned P, typename Tp, unsigned N>
inline const Tp& get(const Ring<Tp, N>& r) { return *std::get<P>(r.addresses()); }

template<unsigned P, typename Tp, unsigned N>
inline Tp& get(Ring<Tp, N>& r) { return *std::get<P>(r.addresses()); }

template<unsigned P, typename Tp, unsigned N>
inline typename Ring<Tp, N>::template ConstAnchor<P> anchor(const Ring<Tp, N>& r)
{ return typename Ring<Tp, N>::template ConstAnchor<P>(r); }

template<unsigned P, typename Tp, unsigned N>
inline typename Ring<Tp, N>::template Anchor<P> anchor(Ring<Tp, N>& r)
{ return typename Ring<Tp, N>::template Anchor<P>(r); }

#endif // SYLVAIN__CODINGAME_INCLUDED