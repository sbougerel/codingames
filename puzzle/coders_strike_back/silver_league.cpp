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
constexpr inline int iabs(int a) {
  int t = a >> (sizeof(int) * 8 - 1);
  return (a^t) - t;
}

// irel() acts as linear rectifier; it returns identity for positive integers,
// but returns 0 for negative integers, such as:
//
//     if (x > 0) return x;
//     else return 0;
constexpr inline int irel(int a) {
  int t = a >> (sizeof(int) * 8 - 1);
  return (a|t) - t;
}

// nirel() returns identity for negative integers, but returns 0 for positive
// integers:
//
//     if (x < 0) return x;
//     else return 0;
constexpr inline int nirel(int a) {
  int t = a >> (sizeof(int) * 8 - 1);
  return (a & t);
}

// isgn() returns `boost` if `gate` is positive or zero, and `-boost` if `gate`
// is negative for integers:
//
//     if (gate >= 0) return boost;
//     if (gate < 0) return -boost;
constexpr inline int isgn(int gate, int boost) {
  int t = gate >> (sizeof(int) * 8 - 1);
  return (boost|t) - t + (t & (-boost));
}

// isgv() returns `pos` if `gate` is positive or zero, and `neg` if `gate` is
// strictly negative for integers:
//
//     if (gate >= 0) return pos;
//     if (gate < 0) return neg;
constexpr inline int isgv(int gate, int pos, int neg) {
  int t = gate >> (sizeof(int) * 8 - 1);
  return (pos|t) - t + (t & neg);
}

// max() returns `a` if `a` is greater than `b`, and `b` otherwise.
constexpr inline int imax(int a, int b) {
  return isgv(a - b, a, b);
}

// min() returns `a` if `a` is lower than `b`, and `b` otherwise.
constexpr inline int imin(int a, int b) {
  return isgv(b - a, a, b);
}

// amp() returns `boost` if `gate` is a positive integer or 0, and returns
// `boost` otherwise for integers in complement 2 notation;
//
//     if (gate >= 0) return boost;
//     else return 0;
constexpr inline int amp(int gate, int boost) {
  int t = gate >> (sizeof(int) * 8 - 1);
  return (boost|t) - t;
}

// namp() returns `boost` if `gate` is a striclty negative integer, and 0
// otherwise for integers in complement 2 notation;
//
//     if (x < 0) return boost;
//     else return 0;
constexpr inline int namp(int gate, int boost) {
  int t = gate >> (sizeof(int) * 8 - 1);
  return (t & boost);
}

// Squares a value.
template<typename Tp>
constexpr inline Tp sq(Tp x) { return x * x; }

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
  int a = iabs(adjacent), o = iabs(opposite);
  int x = a + o;
  x = (sq(x) + S) / (2 * x + 1);
  x = (sq(x) + S) / (2 * x + 1);
  x = (sq(x) + S) / (2 * x + 1);
  return imax(imax(x, a), o);
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
// than 2 degree or error.
//
// Uses a rational polynomial approximation that can be found everywhere on
// Internet, only adapted to integers and degrees. The input range is only
// defined between [0, 1].
//
//     acos(x) ≈ 90 + (ax + bx³) / (1 + cx² + dx⁴)
//
// with:
//     a = -53.807358428
//     b =  52.814341583
//     c = -1.284590624
//     d =  0.295624145
//
// We slightly modify the function above as it returns with [-179, 179] due to
// integer approximations. We instead propose to use:
//
//     acos(x) ≈ 90 + (ax + bx³) / (0.999999f + cx² + dx⁴)
//
// Note that if `hypot` is 0, the result is undefined. `hypot` is assumed to be
// properly computed, and in particular, always greater than `x`.
inline int iacos3(int x, int y, int hypot) {
  constexpr const float A = -53.807358428;
  constexpr const float B = 52.814341583;
  constexpr const float C = -1.284590624;
  constexpr const float D = 0.295624145;
  float f = float(x) / float(hypot);
  float f2 = f * f;
  float q = f * (A + B * f2);
  float d = 0.999999f + f2 * (C + D * f2);
  int r = 90 + int(q / d);
  return isgn(y, r);
}

// Return an angle in degree for the value of the adjacent length `x`, the
// opposite length`y` with less than 0.5% error.
//
inline int iacos2(int x, int y) {
  return iacos3(x, y, ihyp(x, y));
}

struct Vec2 { int x, y; };

constexpr inline int x(const Vec2& a) { return a.x; }
constexpr inline int& x(Vec2& a) { return a.x; }

constexpr inline int y(const Vec2& a) { return a.y; }
constexpr inline int& y(Vec2& a) { return a.y; }

inline std::ostream& operator<< (std::ostream& o, const Vec2& a) {
  return o << "Vec2({" << x(a) << ", " << y(a) << "})";
}

constexpr inline Vec2 operator+ (const Vec2& a, const Vec2& b)
{ return {x(a) + x(b), y(a) + y(b)}; }

constexpr inline Vec2 operator- (const Vec2& a, const Vec2& b)
{ return {x(a) - x(b), y(a) - y(b)}; }

constexpr inline Vec2 operator- (const Vec2& a)
{ return {-x(a), -y(a)}; }

constexpr inline Vec2 operator* (const Vec2& a, int factor)
{ return {x(a) * factor, y(a) * factor}; }

constexpr inline Vec2 operator/ (const Vec2& a, int factor)
{ return {x(a) / factor, y(a) / factor}; }

constexpr inline Vec2 operator* (const Vec2& a, float factor)
{ return {int(float(x(a)) * factor), int(float(y(a)) * factor)}; }

constexpr inline Vec2 operator/ (const Vec2& a, float factor)
{ return {int(float(x(a)) / factor), int(float(y(a)) / factor)}; }

inline Vec2 operator<< (const Vec2& a, int factor)
{ return {x(a) << factor, y(a) << factor}; }

inline Vec2 operator>> (const Vec2& a, int factor)
{ return {x(a) >> factor, y(a) >> factor}; }

constexpr inline bool operator== (const Vec2& a, const Vec2& b)
{ return (x(a) == x(b) && y(a) == y(b)); }

constexpr inline bool operator!= (const Vec2& a, const Vec2& b) { return !(a == b); }

constexpr inline int magsq(const Vec2& a)
{ return sq(x(a)) + sq(y(a)); }

inline int distsq(const Vec2& a, const Vec2& b)
{ return magsq(a - b); }

inline int mag(const Vec2& a) { return ihyp(x(a), y(a)); }

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

constexpr inline const Vec2& low(const Box2& b) { return b.low; }
constexpr inline const Vec2& high(const Box2& b) { return b.high; }
constexpr inline Vec2& low(Box2& b) { return b.low; }
constexpr inline Vec2& high(Box2& b) { return b.high; }

constexpr inline bool within(const Box2& b, const Vec2& v) {
  return (x(v) >= x(low(b)) && x(v) < x(high(b))
          && y(v) >= y(low(b)) && y(v) < y(high(b)));
}

// Given two positions at discreet time t0 and t1 and a distance of closest
// appraoch `sqrad`, perform recursive halving of the time interval to check
// whether the particle collided.
inline int linear_collide(Vec2 x0, Vec2 y0, Vec2 x1, Vec2 y1, int sqrad) {
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

// The basic polar coordinate information, with an angle in degree and a radius
// equivalent to the distance to the pole. Radius is always positive. If it is
// found to be negative, program behaviour will be undefined.
struct Ray2 { int angle, rad; };

constexpr inline int angle(const Ray2& a) { return a.angle; }
constexpr inline int& angle(Ray2& a) { return a.angle; }

constexpr inline int rad(const Ray2& a) { return a.rad; }
constexpr inline int& rad(Ray2& a) { return a.rad; }

constexpr inline Ray2 operator+ (const Ray2& a, const Ray2& b)
{ return {angle(a) + angle(b), rad(a) + rad(b)}; }

constexpr inline Ray2 operator- (const Ray2& a, const Ray2& b)
{ return {angle(a) - angle(b), iabs(rad(a) - rad(b))}; }

constexpr inline Ray2 operator- (const Ray2& a)
{ return {-angle(a), rad(a)}; }

constexpr inline bool operator== (const Ray2& a, const Ray2& b)
{ return (angle(a) == angle(b) && rad(a) == rad(b)); }

constexpr inline bool operator!= (const Ray2& a, const Ray2& b) { return !(a == b); }

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
  int  orient;
  int  rad;
  float mass;
};

// The 4 pairs of accessor in the particule trait are defined below.  A particle
// as a position, a speed, a radius and a mass.
//
constexpr inline const Vec2& pos(const Particle& a) { return a.pos; }
constexpr inline Vec2& pos(Particle& a) { return a.pos; }
constexpr inline const Vec2& spd(const Particle& a) { return a.spd; }
constexpr inline Vec2& spd(Particle& a) { return a.spd; }
constexpr inline int orient(const Particle& a) { return a.orient; }
constexpr inline int& orient(Particle& a) { return a.orient; }
constexpr inline int rad(const Particle& a) { return a.rad; }
constexpr inline int& rad(Particle& a) { return a.rad; }
constexpr inline float mass(const Particle& a) { return a.mass; }
constexpr inline float& mass(Particle& a) { return a.mass; }

constexpr inline bool operator==(const Particle& a, const Particle& b) {
  return (pos(a) == pos(b)
          && spd(a) == spd(b)
          && rad(a) == rad(b)
          && mass(a) == mass(b));
}
constexpr inline bool operator!=(const Particle& a, const Particle& b)
{ return !(a == b); }

inline std::ostream& operator<<(std::ostream& o, const Particle& a) {
  o << "Particle({" << pos(a) << ", " << spd(a) << ", " << orient(a) << ", " << rad(a) << ", " << mass(a) << "})";
  return o;
}

inline Particle linear_motion(const Particle& p) {
  return {spd(p) + pos(p), spd(p), orient(p), rad(p), mass(p)};
}

inline Particle reaction(const Particle& p, const Vec2& t) {
  Vec2 a_ = t / mass(p);
  Vec2 p_ = a_ / 2 + spd(p) + pos(p);
  Vec2 s_ = a_+ spd(p);
  return {p_, s_, orient(p), rad(p), mass(p)};
}

inline Particle reaction(const Particle& p, const Vec2& t, int iterations) {
  Vec2 a_ = t / mass(p);
  Vec2 p_ = (a_ * sq(iterations)) / 2 + spd(p) * iterations + pos(p);
  Vec2 s_ = a_ * iterations + spd(p);
  return {p_, s_, orient(p), rad(p), mass(p)};
}

// ThrustModels functors apply a force on a present particle, computing its
// future position based on the force applied.
//
// InstantThrustModel applies the force as if it the particle had no mass, thus
// propelling it to the expected speed in an instant. This model is not
// realistic, but seem to occur in the puzzles. Mass is ignored in this model.
struct InstantThrustModel
{
  Particle operator() (const Particle& p, const Vec2& t) const {
    Vec2 s_ = t + spd(p);
    Vec2 p_ = s_ + pos(p);
    return {p_, s_, orient(p), rad(p), mass(p)};
  }
};

// RealisticThrustModel applies the force as if it had pushed the particle
// progressively with that force, propelling it to the expected speed in that
// instant. This model is closer to reality, hence the name.
struct RealisticThrustModel
{
  Particle operator() (const Particle& p, const Vec2& t) const {
    return reaction(p, t);
  }
};

// Dragmodels functors return the force of drag execrted on a particle in a
// medium, given its present characteristics.
//
// VaccumDragModel is special in a sense that the model always returns {0, 0}
// for the drag. Good for testing.
struct VaccumDragModel
{
  Vec2 operator() (const Particle&) const { return {0, 0}; }
};

// This simple model seem to be present in several puzzles. For a given type of
// vehicle, it models drag as a ratio of the terminal velocity of the vehicle in
// the medium, applied to the maximum thrust of the vehicle. This ensures that
// the vehicle cannot normaly accelerate beyond its top speed.
template<int MAX_THRUST, int MAX_VELOCITY>
struct BasicDragModel
{
  Vec2 operator() (const Particle& p) const {
    return norm(-spd(p), (mag(spd(p)) * MAX_THRUST) / MAX_VELOCITY);
  }
};

// CoastingAction just let the particle decelrate by drag. Important to compute
// break distance under drag in any phyical model.
struct CoastingAction
{
  Ray2 operator() (const Particle& p) const { return {orient(p), 0}; }
};

// ConstantAction just makes the particle accelerate with a constant thrust
// applied in the same direction. Good for tests.
struct ConstantAction
{
  ConstantAction(const Vec2& thrust) : _thrust(ray(thrust)) { }
  Ray2 operator() (const Particle&) const { return _thrust; }
private:
  Ray2 _thrust;
};

// TargetAction just makes the particle move toward a target with a constant
// thrust. It is only slightly more useful than BasicModel in puzzles. The
// target is set in the constructor.
struct TargetAction
{
  TargetAction(const Vec2& target, int thrust)
    : _target(target), _thrust(thrust) { }
  Ray2 operator() (const Particle& p) const {
    return {angle(ray(_target - pos(p))), _thrust};
  }
private:
  Vec2 _target;
  int _thrust;
};

// AdvTargetAction acts like TargetAction but also tries to compensate its own
// lateral motion to reach the target in a straighter line. This model is only a
// few lines but can be used to simulate advanced racing bots.
template<int MAX_THRUST, int MAX_CORRECTION> // maximum correction angle
struct AdvTargetAction
{
  AdvTargetAction(const Vec2& target, int radius) : _target(target), _radius(radius) { }
  Ray2 operator() (const Particle& p) const {
    constexpr const int MAX_COMP_ANGLE = 90;  // compensate when facing
    constexpr const int ACCEL_ANGLE    = 100; // start accelerating
    if (magsq(spd(p)) < 100)
      return TargetAction(_target, MAX_THRUST)(p);
    Ray2 pro    = ray(spd(p));
    if (magsq(spd(p) + pos(p) - _target) < sq(_radius)
        && iabs(angle(pro) - orient(p)) < MAX_CORRECTION)
      { return {angle(pro), MAX_THRUST}; }
    Ray2 dir    = ray(_target - pos(p));
    int  pro_d  = angle(dir) - angle(pro);
    Ray2 push   = dir;
    if (iabs(pro_d) < MAX_COMP_ANGLE)
      { angle(push) = angle(dir) + isgn(pro_d, imin(iabs(pro_d), MAX_CORRECTION)); }
    int  ori_d  = iabs(angle(push) - orient(p));
    if (ori_d > ACCEL_ANGLE) { rad(push) = 0; }
    else if (ori_d < MAX_CORRECTION * 2) { rad(push) = MAX_THRUST; }
    else { rad(push) = isin(((ACCEL_ANGLE - ori_d) * 90) / (ACCEL_ANGLE - MAX_CORRECTION * 2), MAX_THRUST + 1); }
    return push;
  }
private:
  Vec2 _target;
  int _radius;
};

// Physics are modeled with a ThrustModel and a DragModel.
template<typename ThrustModel, typename DragModel>
struct Physics : private ThrustModel, DragModel {
  Physics(const ThrustModel& tm = ThrustModel(),
          const DragModel& dm = DragModel())
    : ThrustModel(tm), DragModel(dm) { }
  const ThrustModel& thrustModel() const { return *this; }
  const DragModel& dragModel() const { return *this; }
};

// `reaction`, `iterate_reaction` and `until_reaction` project actions on
// particles to compute the future of a particle based on its
// known present and a phyical model.
template<typename ThrustModel, typename DragModel>
inline Particle reaction(const Particle& p, const Vec2& t,
                         const Physics<ThrustModel, DragModel>& phy) {
  return phy.thrustModel()(p, t + phy.dragModel()(p));
}

template<typename Action, typename ThrustModel, typename DragModel>
inline Particle iterate_reaction(unsigned times, Particle p, const Action& a,
                                 const Physics<ThrustModel, DragModel>& phy) {
  for (unsigned i = 0; i < times; ++i) { p = reaction(p, vec(a(p)), phy); }
  return p;
}

template<typename Action, typename ThrustModel, typename DragModel, typename Predicate>
inline Particle until_reaction(Particle p, const Action& a, const Predicate& t,
                               const Physics<ThrustModel, DragModel>& phy) {
  while (!t(p)) { p = reaction(p, vec(a(p)), phy); }
  return p;
}

// Collision Detection algorithm. Returns the distance of collision, an a
// posteriori estimate of the closest approach between 2 particles (squared)
// and the time to collision. When closest approach <= distance of collision,
// collision has occured.
//
// The estimation stops after max_iter, which is 100 by default, or when one of
// the particle gets out of the bounding box.
//
// The estimation is based on the movement Model associated with both
// particles. At each turn, the model is updated with the particles' positions,
// and it queries the thrust for each of the particles.
//
template<typename ActionA, typename ActionB, typename Physics>
inline std::tuple<int, int, int>
collide_two(Particle a0, Particle b0, const ActionA& ma, const ActionB& mb,
            const Physics& phy, int max_iter = 100,
            const Box2& bb = {{-10000,-10000}, {10000, 10000}}) {
  int sqrad = sq(rad(a0)) + sq(rad(b0));
  int best_approach = distsq(pos(a0), pos(b0));
  if (best_approach <= sqrad)
    { return std::make_tuple(sqrad, best_approach, 0); }
  int i = 0;
  for (; i < max_iter; ++i) {
    Particle a1 = reaction(a0, ma, phy);
    Particle b1 = reaction(b0, mb, phy);
    if (!within(bb, pos(a1)) || !within(bb, pos(b1))) break;
    int approach = linear_collide(pos(a0), pos(b0), pos(a1), pos(b1), sqrad);
    if (approach <= sqrad)
      { return std::make_tuple(sqrad, approach, i); }
    if (approach < best_approach) { best_approach = approach; }
    a0 = a1; b0 = b1;
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
    Tp* copy = _addr[N - 1];
    for (unsigned i = N - 1; i > 0; --i) { _addr[i] = _addr[i - 1]; }
    _addr[0] = copy;
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

#include <vector>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::tuple;
using std::vector;

constexpr const int MAX_THRUST       = 100;
constexpr const int MAX_SPEED        = 660; // fast compute drag
constexpr const int POD_RADIUS       = 400;
constexpr const int POD_MASS         = 1;
constexpr const int CP_RADIUS        = 600;
constexpr const int MAP_WIDTH        = 16000;
constexpr const int MAP_HEIGHT       = 9000;
constexpr const int MAP_SEMI_WIDTH   = MAP_WIDTH / 2;
constexpr const int MAP_SEMI_HEIGHT  = MAP_HEIGHT / 2;
constexpr const int MAX_POD_ROTATION = 18; // rotation speed in degree
const char BOOST[]                   = "BOOST";
const char SHIELD[]                  = "SHIELD";

constexpr inline Vec2 to_centered(const Vec2& v) {
  return v + Vec2{-MAP_SEMI_WIDTH, -MAP_SEMI_HEIGHT};
}

constexpr inline Vec2 to_local(const Vec2& v) {
  return v + Vec2{MAP_SEMI_WIDTH, MAP_SEMI_HEIGHT};
}
struct State {
    Particle myPod;
    Particle thPod;
    Vec2     myCpPos;
    Ray2     myCpRay;
};

inline State readState() {
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
  return State{Particle{to_centered({x, y}), {0, 0}, 0, POD_RADIUS, POD_MASS},
               Particle{to_centered({opponentY, opponentY}), {0, 0}, 0, POD_RADIUS, POD_MASS},
               to_centered({nextCheckpointX, nextCheckpointY}),
               Ray2{nextCheckpointAngle, nextCheckpointDist}};
}

inline void updateState(State& curr, const State& prev) {
  spd(curr.myPod) = pos(curr.myPod) - pos(prev.myPod);
  spd(curr.thPod) = pos(curr.thPod) - pos(prev.thPod);
  // Unit circle goes clockwise, while orientation is anti-clockwise
  orient(curr.myPod) = angle(ray(curr.myCpPos - pos(curr.myPod))) - angle(curr.myCpRay);
}

typedef vector<tuple<Vec2, bool>> CheckPoints;

typedef Ring<State, 3> History;

inline void thrust(int x, int y, int t) {
  cout << x << " " << y << " " << imin(t, MAX_THRUST) << endl;
}
inline void thrust(const Vec2& p, int t) {
  Vec2 l = to_local(p);
  thrust(x(l), y(l), t);
}

inline void boost(int x, int y) {
    cout << x << " " << y << " BOOST" << endl;
}
inline void boost(const Vec2& p) {
  Vec2 l = to_local(p);
  boost(x(l), y(l));
}

inline void shield(int x, int y) {
    cout << x << " " << y << " SHIELD" << endl;
}
inline void shield(const Vec2& p) {
  Vec2 l = to_local(p);
  shield(x(l), y(l));
}

/**
 * Rotation & acceleration test
 **/
int main()
{
  Physics<InstantThrustModel, BasicDragModel<MAX_THRUST, MAX_SPEED>> phys;
  bool boost_used = false;
  History hist(State{Particle{{0, 0}, {0, 0}, 0, POD_RADIUS, POD_MASS},
                     Particle{{0, 0}, {0, 0}, 0, POD_RADIUS, POD_MASS},
                     {0, 0}, {0, 0}});
  auto curr = anchor<0>(hist);
  auto prev = anchor<1>(hist);
  *curr = readState();
  Ray2 push = AdvTargetAction<MAX_THRUST, MAX_POD_ROTATION>(curr->myCpPos, CP_RADIUS - 50)(curr->myPod);
  if (linear_collide(pos(curr->myPod), pos(curr->thPod),
                     pos(reaction(curr->myPod, vec(push), phys)),
                     pos(curr->thPod) + spd(curr->thPod),
                     sq(POD_RADIUS * 2)) <= sq(POD_RADIUS * 2))
    {
      shield(pos(curr->myPod) + vec({angle(push), 2000}));
    }
  if (iabs(angle(ray(spd(curr->myPod))) - angle(ray(curr->myCpPos - pos(curr->myPod)))) < MAX_POD_ROTATION
      && iabs(angle(curr->myCpRay)) < MAX_POD_ROTATION
      && rad(curr->myCpRay) > 2000
      && rad(push) == MAX_THRUST) {
    boost(pos(curr->myPod) + vec({angle(push), 2000}));
    boost_used = true;
  }
  else
    thrust(pos(curr->myPod) + vec({angle(push), 2000}), rad(push));

  // game loop
  while (1) {
    hist.rotate();
    *curr = readState();
    updateState(*curr, *prev);
    Ray2 push = AdvTargetAction<MAX_THRUST, MAX_POD_ROTATION>(curr->myCpPos, CP_RADIUS - 50)(curr->myPod);
    if (linear_collide(pos(curr->myPod), pos(curr->thPod),
                       pos(reaction(curr->myPod, vec(push), phys)),
                       pos(curr->thPod) + spd(curr->thPod),
                       sq(POD_RADIUS * 2)) <= sq(POD_RADIUS * 2))
      {
        shield(pos(curr->myPod) + vec({angle(push), 2000}));
      }
    if (!boost_used
        && iabs(angle(ray(spd(curr->myPod))) - angle(ray(curr->myCpPos - pos(curr->myPod)))) < MAX_POD_ROTATION
        && iabs(angle(curr->myCpRay)) < MAX_POD_ROTATION
        && rad(curr->myCpRay) > 2000
        && rad(push) == MAX_THRUST) {
      boost(pos(curr->myPod) + vec({angle(push), 2000}));
      boost_used = true;
    }
    else
      thrust(pos(curr->myPod) + vec({angle(push), 2000}), rad(push));
    cerr << "Last pos " << pos(prev->myPod) << " Curr pos " << pos(curr->myPod) << endl;
    cerr << "Speed " << spd(curr->myPod) << " (" << mag(spd(curr->myPod)) << ")" << endl;
    cerr << "Accel " << spd(curr->myPod) - spd(prev->myPod) << " (" << mag(spd(curr->myPod) - spd(prev->myPod)) << ")" << endl;
    cerr << "Rotated " << angle(curr->myCpRay) - angle(prev->myCpRay) << endl;
  }
}
