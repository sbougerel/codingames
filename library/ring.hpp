// Ring & Anchor are 2 simple objects that store objects in a contiguous location
// and then rotate addressing to each objects stored.
//
// Ring is designed to store a small amount of large types. It is non-copyable
// but movable.
#ifndef SYLVAIN__CODINGAME_RING
#define SYLVAIN__CODINGAME_RING

#include <array>
#include <vector>

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


#endif // SYLVAIN__CODINGAME_RING
