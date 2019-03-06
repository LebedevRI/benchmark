#include "benchmark/benchmark.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

// Future preliminary libc++ code, from Marshall Clow.
namespace std {
template <class _Tp>
__inline _Tp midpoint(_Tp __a, _Tp __b) noexcept {
  using _Up = typename std::make_unsigned<typename remove_cv<_Tp>::type>::type;

  int __sign = 1;
  _Up __m = __a;
  _Up __M = __b;
  if (__a > __b) {
    __sign = -1;
    __m = __b;
    __M = __a;
  }
  return __a + __sign * _Tp(_Up(__M - __m) >> 1);
}
}  // namespace std

template <typename T>
std::vector<T> getVectorOfRandomNumbers(size_t count) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<T> dis(std::numeric_limits<T>::min(),
                                       std::numeric_limits<T>::max());
  std::vector<T> v;
  v.reserve(count);
  std::generate_n(std::back_inserter(v), count,
                  [&dis, &gen]() { return dis(gen); });
  assert(v.size() == count);
  return v;
}

struct RandRand {
  template <typename T>
  static std::pair<std::vector<T>, std::vector<T>> Gen(size_t count) {
    return std::make_pair(getVectorOfRandomNumbers<T>(count),
                          getVectorOfRandomNumbers<T>(count));
  }
};
struct ZeroRand {
  template <typename T>
  static std::pair<std::vector<T>, std::vector<T>> Gen(size_t count) {
    return std::make_pair(std::vector<T>(count, T(0)),
                          getVectorOfRandomNumbers<T>(count));
  }
};

template <class T, class Gen>
void BM_StdMidpoint(benchmark::State& state) {
  const size_t Length = state.range(0);

  const std::pair<std::vector<T>, std::vector<T>> Data =
      Gen::template Gen<T>(Length);
  const std::vector<T>& a = Data.first;
  const std::vector<T>& b = Data.second;
  assert(a.size() == Length && b.size() == a.size());

  benchmark::ClobberMemory();
  benchmark::DoNotOptimize(a);
  benchmark::DoNotOptimize(a.data());
  benchmark::DoNotOptimize(b);
  benchmark::DoNotOptimize(b.data());

  for (auto _ : state) {
    for (size_t i = 0; i < Length; i++) {
      const auto calculated = std::midpoint(a[i], b[i]);
      benchmark::DoNotOptimize(calculated);
    }
  }
  state.SetComplexityN(Length);
  state.counters["midpoints"] =
      benchmark::Counter(Length, benchmark::Counter::kIsIterationInvariant);
  state.counters["midpoints/sec"] =
      benchmark::Counter(Length, benchmark::Counter::kIsIterationInvariantRate);
  const size_t BytesRead = 2 * sizeof(T) * Length;
  state.counters["bytes_read/iteration"] =
      benchmark::Counter(BytesRead, benchmark::Counter::kDefaults,
                         benchmark::Counter::OneK::kIs1024);
  state.counters["bytes_read/sec"] = benchmark::Counter(
      BytesRead, benchmark::Counter::kIsIterationInvariantRate,
      benchmark::Counter::OneK::kIs1024);
}

template <typename T>
static void CustomArguments(benchmark::internal::Benchmark* b) {
  const size_t L2SizeBytes = 2 * 1024 * 1024;
  // What is the largest range we can check to always fit within given L2 cache?
  const size_t MaxLen = L2SizeBytes / /*total bufs*/ 2 /
                        /*maximal elt size*/ sizeof(T) / /*safety margin*/ 2;
  b->RangeMultiplier(2)->Range(1, MaxLen)->Complexity(benchmark::oN);
}

// Both of the values are random.
// The comparison is unpredictable.
BENCHMARK_TEMPLATE(BM_StdMidpoint, int32_t, RandRand)
    ->Apply(CustomArguments<int32_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, uint32_t, RandRand)
    ->Apply(CustomArguments<uint32_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, int64_t, RandRand)
    ->Apply(CustomArguments<int64_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, uint64_t, RandRand)
    ->Apply(CustomArguments<uint64_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, int16_t, RandRand)
    ->Apply(CustomArguments<int16_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, uint16_t, RandRand)
    ->Apply(CustomArguments<uint16_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, int8_t, RandRand)
    ->Apply(CustomArguments<int8_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, uint8_t, RandRand)
    ->Apply(CustomArguments<uint8_t>);

// One value is always zero, and another is bigger or equal than zero.
// The comparison is predictable.
BENCHMARK_TEMPLATE(BM_StdMidpoint, uint32_t, ZeroRand)
    ->Apply(CustomArguments<uint32_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, uint64_t, ZeroRand)
    ->Apply(CustomArguments<uint64_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, uint16_t, ZeroRand)
    ->Apply(CustomArguments<uint16_t>);
BENCHMARK_TEMPLATE(BM_StdMidpoint, uint8_t, ZeroRand)
    ->Apply(CustomArguments<uint8_t>);
