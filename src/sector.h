#pragma once

#include <iostream>
#include <vector>
#include <numeric>

#include "yaml-cpp/yaml.h"

#include "utils.h"


class Reduce;

class RawIntegral {
public:
  RawIntegral() = default;

  explicit RawIntegral(unsigned n) : indices(n) {}

  explicit RawIntegral(unsigned n, int i) : indices(n, i) {}

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

  friend bool operator==(const RawIntegral &lhs, const RawIntegral &rhs) {
    return lhs.indices == rhs.indices;
  }

  friend bool operator<(const RawIntegral &lhs, const RawIntegral &rhs) {
    return lhs.indices < rhs.indices;
  }

  friend std::ostream &operator<<(std::ostream &os, const RawIntegral &integral) {
    os << integral.indices;
    return os;
  }

private:
  std::vector<int> indices;
};

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

class Sector {
public:
  friend class Reduce;

private:
  // sector number
  unsigned _id = 0;
  // number of targets
  unsigned _ntargets = 0;
  // maximum depth
  unsigned _depth = 0;
  // maximum rank
  unsigned _rank = 0;
  // super sectors
  std::vector<unsigned> _superSectors;
  // sub sectors
  std::vector<unsigned> _subSectors;
};