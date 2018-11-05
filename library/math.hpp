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
  inline int do_sin(int angle, int scale) {
    constexpr const int FACTOR_1 = 57; // similar to 180 / M_PI
    constexpr const int FACTOR_2 = 6 * sq<FACTOR_1>();
    return (((scale * angle) / FACTOR_1) * (FACTOR_2 - sq(angle))) / (FACTOR_2);
  }

  inline int do_cos(int angle, int scale) {
    constexpr const int FACTOR_1 = 57; // similar to 180 / M_PI
    constexpr const int FACTOR_2 = 2 * sq<FACTOR_1>();
    return isgn(angle, (scale * FACTOR_2 - scale * sq(90 - abs(angle))) / FACTOR_2);
  }
}

// Return the sine value for an `angle` in degree, multipled by arbitrary
// precision. This method is only defined for the interval within (-450, 450).
//
// It uses only the first degree of the Taylor serie for either sin or cos,
// since the first degree is only precise until (-45, 45). It should have less
// that 0.5% of error.
//
//   1. map (-450, 450) to (-90, 90) with rectifiers
//   2. use cos or sine first taylor order
inline int isin(int angle, int scale) {
  angle = angle + 2 * (- irel(angle - 90) - nirel(angle + 90)
                       + irel(angle - 270) + nirel(angle + 270));
  return isgv(abs(angle) - 45, details::do_cos(angle, scale), details::do_sin(angle, scale));
}

// This method is only defined for the interval within (-360, 360).
inline int icos(int angle, int scale) {
  return isin(90 - angle, scale);
}
