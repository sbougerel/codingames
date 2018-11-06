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

inline int sq(int x) { return x * x; }

template<int N>
constexpr int sq() { return N*N; }

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

namespace details {
  template<int Factor>
  inline int do_sin(int angle, int scale) {
    constexpr const int Factor2 = 6 * sq<Factor>();
    return (((scale * angle) / Factor) * (Factor2 - sq(angle))) / Factor2;
  }

  template<int Factor>
  inline int do_cos(int angle, int scale) {
    constexpr const int Factor2 = 2 * sq<Factor>();
    return (scale * Factor2 - scale * sq(angle)) / Factor2;
  }
}

// Return the sine value for an `angle` in degree, multipled by arbitrary
// precision. This method is only defined for the interval within (-450, 450).
//
// It uses only the first degree of the Taylor serie for either sin or cos,
// since the first degree is only precise until (-45, 45). It should have less
// that 0.5% of error.
//
inline int isin_(int angle, int scale) {
  angle = angle + 2 * (- irel(angle - 90) - nirel(angle + 90)
                       + irel(angle - 270) + nirel(angle + 270));
  // 57 ~= 180 / M_PI
  return isgv(abs(angle) - 45, isgn(angle, details::do_cos<57>(90 - abs(angle), scale)), details::do_sin<57>(angle, scale));
}

// Return the sine value for an `angle` in degree, multipled by arbitrary
// precision.
//
// It does it by mapping degree angles into into 512 values. The first 7 bits
// are the angle within a 45 degree quadrant. The next 3 bits control the choice
// of sign and operations to perform: sin(x) or cos(x) are each more precise
// when x is close to 0.
//
// It should have less than 1% of error.
inline int isin(int angle, int scale) {
  int aa = ( abs(angle) * 128 + 45 ) / 90;
  int la = aa & 0x3F; // angle part
  int h = ((aa & 0x100) >> 8) - 1; // half Q1-4 or Q5-8: 0 or -1
  int op = (aa & 0xA0) >> 6;
  // 81 = 256 / M_PI
  int r;
  switch (op) {
  case 0: r = details::do_sin<81>(la, scale); break;
  case 1: r = details::do_cos<81>(64 - la, scale); break;
  case 2: r = details::do_cos<81>(la, scale); break;
  case 3: r = details::do_sin<81>(64 - la, scale); break;
  }
  return isgn(angle, h * r);
}

// This method is only defined for the interval within (-360, 360).
inline int icos(int angle, int scale) {
  return isin(90 - angle, scale);
}
