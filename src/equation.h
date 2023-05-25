#pragma once

#include <map>
#include <algorithm>

#include "arith/umod.h"

// equation over finite field
class EquationFF {
public:
  umod64 operator[](unsigned i) const {
    return _eq.at(i);
  }

  umod64 &operator[](unsigned i) {
    return _eq[i];
  }

  // insert a new item
  void insert(unsigned integral, umod64 coeff) {
    _eq.insert({integral, coeff});
  }

  // number of items
  [[nodiscard]] unsigned size() const {
    return _eq.size();
  }

  // get the first integral number
  unsigned first_integral() {
    return _eq.begin()->first;
  }

  // get the first coefficient
  umod64 first_coeff() {
    return _eq.begin()->second;
  }

  // is empty
  [[nodiscard]] bool empty() const {
    return _eq.empty();
  }

  // clear zero items
  void erase_zero() {
    std::erase_if(_eq, [](const auto &item) { return item.second == 0; });
  }

  // normalize: set the first coeff to one
  void normalize() {
    umod64 scale = _eq.begin()->second;
    for (auto &item: _eq)
      item.second /= scale;
  }

  // gauss elimination: eliminate the other equation from this one
  // it is assumed that the other equation has been normalized
  void eliminate(const EquationFF &other) {
    umod64 scale = _eq.begin()->second;
    for (const auto &item: other._eq)
      _eq[item.first] -= item.second * scale;
    this->erase_zero();
  }

  // get the underline eq, only readable
  const auto &eq() {
    return _eq;
  }

private:
  std::map<unsigned, umod64, std::greater<>> _eq;
};
