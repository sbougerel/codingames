// It wouldn't be fast integer math if we didn't guarantee that things run
// faster that normal double operations.

#define BOOST_TEST_MODULE Approximate Integer Performance
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <chrono>
#include <iostream>

#include "included.hpp"

template<typename Tp>
constexpr inline void silence(Tp) {}

struct random_int {
  int seed() {
    std::srand(std::time(nullptr));
    return 0;
  }
  int operator() () {
    static int unused = seed();
    silence(unused);
    return std::rand();
  }
};

template<size_t N>
inline std::vector<int> random_ints() {
  std::vector<int> v(N);
  random_int r;
  for (size_t i = 0; i < N; ++i) v[i] = r();
  return v;
}

BOOST_AUTO_TEST_CASE(test_abs){
  constexpr const int N = 10000000;

  // std::abs() is actually really fast. Here we just test my version is not too
  // slow since sgn(), rectify(), etc. all use the same technique.
  std::vector<int> shuffle(random_ints<N>());

  std::chrono::duration<double> elapsed_std;
  {
    int j = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i : shuffle) { j += std::abs(i); }
    auto end = std::chrono::high_resolution_clock::now();
    elapsed_std = end-start;
    std::cout << "std::abs() elapsed time for " << N << " iterations:\t"
              << elapsed_std.count() << "s (" << j << ")" << std::endl;
  }

  std::chrono::duration<double> elapsed_my;
  {
    int j = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i : shuffle) { j += iabs(i); }
    auto end = std::chrono::high_resolution_clock::now();
    elapsed_my = end-start;
    std::cout << "iabs() elapsed time for " << N << " iterations:\t"
              << elapsed_my.count() << "s (" << j << ")" << std::endl;
  }

  BOOST_CHECK_GT(elapsed_std.count() * 1.5, elapsed_my.count());
}

BOOST_AUTO_TEST_CASE(test_sine){
  constexpr const int N = 10000000;

  // again I find that std::sin is impressively fast, for something so
  // accurate.
  std::vector<int> deg_shuffle(random_ints<N>());
  for (size_t i= 0; i < N; ++i) deg_shuffle[i] = deg_shuffle[i] % 720;
  std::vector<double> rad_shuffle;
  for (size_t i= 0; i < N; ++i) rad_shuffle.push_back(double(deg_shuffle[i]) * M_PI / 180.0);

  std::chrono::duration<double> elapsed_base;
  {
    int j = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i : deg_shuffle) { j += i; }
    auto end = std::chrono::high_resolution_clock::now();
    elapsed_base = end-start;
    std::cout << "base iteration time for " << N << " iterations:\t"
              << elapsed_base.count() << "s (" << j << ")" << std::endl;
  }

  std::chrono::duration<double> elapsed_std;
  {
    double j = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (double i : rad_shuffle) { j += std::sin(i); }
    auto end = std::chrono::high_resolution_clock::now();
    elapsed_std = end-start;
    std::cout << "std::sin() elapsed time for " << N << " iterations:\t"
              << elapsed_std.count() << "s (" << j << ")" << std::endl;
  }

  std::chrono::duration<double> elapsed_my;
  {
    int j = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i : deg_shuffle) { j += isin(i, 1000); }
    auto end = std::chrono::high_resolution_clock::now();
    elapsed_my = end-start;
    std::cout << "isin() elapsed time for " << N << " iterations:\t"
              << elapsed_my.count() << "s (" << j << ")" << std::endl;
  }

  BOOST_CHECK_GT(elapsed_std.count(), elapsed_my.count());
}
