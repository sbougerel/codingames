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

template<int N>
constexpr int sq() { return N*N; }
inline int sq(int x) { return x * x; }

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
  // Curve-fitted 3rd order Taylor approximation to the sine function that only
  // takes 256 positive angles values representing angles between 0 and PI.
  inline int do_sin(int angle, int scale) {
    constexpr const int Factor = 81;                // 256 / PI ~= 81
    constexpr const int Factor2 = 6 * sq<Factor>(); // 3! * (256/PI ~= 81)^2
    return (((scale * angle) / Factor) * (Factor2 - (sq(angle) * 8) / 9)) / Factor2;
  }
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
  int s = angle >> (sizeof(int) * 8 - 1);     // sign mask
  int aa = (((angle^s) - s) * 128 + 45) / 90; // absolute angle
  int la = aa & 0x7F;                         // quadrant angle
  int h = -((aa & 0x100) >> 8);               // first or second half: 0 or -1
  int q = -((aa & 0x80)  >> 7);               // first or second quadrant: 0 or -1
  int r = (details::do_sin(la, scale) | q) + (details::do_sin(128 - la, scale) & q) - q;
  s = h^s;                                    // XOR halves and sign together
  return (r^s) - s;
}

inline int icos(int angle, int scale) {
  return isin(90 - angle, scale);
}
