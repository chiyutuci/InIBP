#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <numeric>

#include "yaml-cpp/yaml.h"
#include "ginac/ginac.h"
#include "arith/umod.h"

#include "utils.h"
#include "equation.h"

class Reduce;

class RawIntegral {
public:
  RawIntegral() = default;

  explicit RawIntegral(unsigned n) : indices(n) {}

  explicit RawIntegral(unsigned n, int i) : indices(n, i) {}

  explicit RawIntegral(const std::vector<bool> &lines) : indices(lines.size()) {
    for (unsigned i = 0; i < lines.size(); ++i)
      indices[i] = lines[i];
  }

  explicit RawIntegral(const std::vector<int> &v) : indices(v) {}

  explicit RawIntegral(std::vector<int> &&v) : indices(std::move(v)) {}

  int operator[](unsigned i) const { return indices[i]; }

  int &operator[](unsigned i) { return indices[i]; }

  // length of propagators
  [[nodiscard]] unsigned size() const { return indices.size(); }

  // sum of positive indices
  [[nodiscard]] unsigned depth() const {
    return std::accumulate(indices.begin(), indices.end(), 0, [](int a, int b) {
      return b > 0 ? a + b : a;
    });
  }

  // sum of negative indices
  [[nodiscard]] unsigned rank() const {
    return std::accumulate(indices.begin(), indices.end(), 0, [](int a, int b) {
      return b < 0 ? a - b : a;
    });
  }

  // sector of the integral
  [[nodiscard]] unsigned sector() const {
    unsigned sector = 0;
    for (unsigned i = 0; i < indices.size(); ++i)
      if (indices[i] > 0)
        sector |= 1 << i;
    return sector;
  }

  bool operator==(const RawIntegral &other) const {
    return other.indices == indices;
  }

  bool operator!=(const RawIntegral &other) const {
    return indices < other.indices;
  }

  friend bool operator<(const RawIntegral &lhs, const RawIntegral &rhs) {
    return lhs.indices < rhs.indices;
  }

  friend RawIntegral operator+(const RawIntegral &lhs, const RawIntegral &rhs) {
    RawIntegral integral(lhs.size());
    for (unsigned i = 0; i < integral.size(); ++i)
      integral[i] = lhs[i] + rhs[i];
    return integral;
  }

  friend RawIntegral operator-(const RawIntegral &lhs, const RawIntegral &rhs) {
    RawIntegral integral(lhs.size());
    for (unsigned i = 0; i < integral.size(); ++i)
      integral[i] = lhs[i] - rhs[i];
    return integral;
  }

  friend std::ostream &operator<<(std::ostream &os, const RawIntegral &integral) {
    os << integral.indices;
    return os;
  }

  friend class std::hash<RawIntegral>;

private:
  std::vector<int> indices;
};

namespace std {
  template<>
  struct hash<RawIntegral> {
    std::size_t operator()(const RawIntegral &integral) const {
      std::size_t hash = integral.size();
      for (auto i: integral.indices)
        hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
      return hash;
    }
  };
}

namespace YAML {
  template<>
  struct convert<RawIntegral> {
    static bool decode(const Node &node, RawIntegral &integral) {
      if (!node.IsSequence())
        return false;
      integral = RawIntegral(node.as<std::vector<int>>());
      return true;
    }
  };
}

// ibp relation
// first:  integral indices
// second: coefficient
typedef std::vector<std::pair<RawIntegral, GiNaC::ex>> IBPProto;
// ibp relation over finite field
// first:  integral indices
// second: coefficients of indices
typedef std::vector<std::pair<RawIntegral, std::vector<umod64>>> IBPProtoFF;

class Sector {
public:
  friend class Reduce;

  // generate seeds and read targets
  void prepare_targets(const std::vector<RawIntegral> &);
  // run the reduction
  void run_reduce(const std::vector<IBPProtoFF> &);

private:
  // generate seeds satisfying the depth and rank
  void _generate_seeds();

  // generate an equation which contains a certain integral
  void _generate_equation(unsigned, const std::vector<IBPProtoFF> &);

public:
  // dp of combinations
  // combinations[<number, sum>] -> [combinations]
  static std::map<std::pair<int, int>, std::vector<std::vector<int>>> combinations;
  // fill the combinations up to <number, sum>
  static void generate_combinations(int number, int sum);

private:
  // sector number
  unsigned _id = 0;
  // number of propagators
  unsigned _nprops = 0;
  // maximum depth
  unsigned _depth = 0;
  // maximum rank
  unsigned _rank = 0;

  // lines
  std::vector<bool> _lines;
  // super sectors
  std::vector<unsigned> _superSectors;
  // sub sectors
  std::vector<unsigned> _subSectors;
  // seeds: weight to integral
  std::vector<RawIntegral> _seeds;
  // seeds: integral to weight
  std::unordered_map<RawIntegral, unsigned> _weights;

  // target integrals
  std::set<unsigned, std::greater<>> _targets;

  // auxiliary targets
  std::queue<unsigned> _auxTargets;
  // used targets
  std::set<unsigned, std::greater<>> _usedTargets;
  // used ibp equations
  // first: integral
  // second: number of ibp
  std::set<std::pair<unsigned, unsigned>> _usedIBP;
  // the ibp system
  std::vector<EquationFF> _systemFF;
  // line number of each pivot
  std::unordered_map<unsigned, unsigned> _lineNumber;
};