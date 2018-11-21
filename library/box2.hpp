#ifndef BOX_HPP
#define BOX_HPP

#include "vec2.hpp"

struct Box2 { Vec2 low, high; };

inline const Vec2& low(const Box2& b) { return b.low; }
inline const Vec2& high(const Box2& b) { return b.high; }
inline Vec2& low(Box2& b) { return b.low; }
inline Vec2& high(Box2& b) { return b.high; }

inline bool within(const Box2& b, const Vec2& v) {
  return (x(v) >= x(low(b)) && x(v) < x(high(b))
          && y(v) >= y(low(b)) && y(v) < y(high(b)));
}

#endif
