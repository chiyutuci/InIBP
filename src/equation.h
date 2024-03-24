#pragma once

#include <map>
#include <algorithm>
#include <iostream>

#include "utils.h"
#include "arith/umod.h"
#include "ginac/ginac.h"

// equation over finite field
class EquationFF {
public:
  unsigned operator[](unsigned i) const {
    return _eq[i].first;
  }

  unsigned &operator[](unsigned i) {
    return _eq[i].first;
  }

 [[nodiscard]] umod64 coeff(unsigned i) const {
    return _eq[i].second;
  }

  // insert a new item
  void insert(unsigned integral, umod64 coeff) {
    _eq.emplace_back(integral, coeff);
  }

  // number of items
  [[nodiscard]] unsigned size() const {
    return _eq.size();
  }

  // get the first integral number
  [[nodiscard]] unsigned first_integral() const {
    return _eq[0].first;
  }

  // get the first coefficient
  umod64 first_coeff() {
    return _eq[0].second;
  }

  // is empty
  [[nodiscard]] bool empty() const {
    return _eq.empty();
  }

  void sort() {
    std::sort(_eq.begin(), _eq.end(),
              [](const std::pair<unsigned, umod64>& i, const std::pair<unsigned, umod64>& j){
                  return i.first > j.first;
    });
  }

  // clear zero items
  void erase_zero() {
    std::erase_if(_eq, [](const auto &item) { return item.second == 0; });
  }

  // normalize: set the first coeff to one
  void normalize() {
    umod64 scale = _eq[0].second;
    for (auto &item: _eq)
      item.second /= scale;
  }

  // gauss elimination: eliminate the other equation from this one
  // it is assumed that the other equation has been normalized
  void eliminate(const EquationFF &other, unsigned index) {
    std::vector<std::pair<unsigned, umod64>> eq;

    umod64 scale = _eq[index].second;
    unsigned iother = 0, ithis = 0;
    while (iother < other.size() && ithis < _eq.size()) {
      if(other[iother] > _eq[ithis].first) {
        eq.emplace_back(other[iother], -scale * other._eq[iother].second);
        ++iother;
      }
      else if (other[iother] == _eq[ithis].first) {
        umod64 coeff = _eq[ithis].second - scale * other.coeff(iother);
        if (coeff != 0)
          eq.emplace_back(_eq[ithis].first, coeff);
        ++iother;
        ++ithis;
      }
      else {
        eq.emplace_back(_eq[ithis].first, _eq[ithis].second);
        ++ithis;
      }
    }

    if (iother < other.size())
      for(; iother < other.size(); ++iother)
        eq.emplace_back(other[iother], -scale * other.coeff(iother));
    if (ithis < _eq.size())
      for (; ithis < _eq.size(); ++ithis)
        eq.emplace_back(_eq[ithis]);

    std::swap(eq, _eq);
  }

  // get the underline eq, only readable
  const auto &eq() {
    return _eq;
  }

  // order of equations
  bool operator<(const EquationFF& other) const {
    if (this->first_integral() != other.first_integral())
      return this->first_integral() < other.first_integral();
    else {
      if (this->size() != other.size())
        return this->size() < other.size();
      else
        return this->eqnum < other.eqnum;
    }
  }

private:
  std::vector<std::pair<unsigned, umod64>> _eq;
public:
  unsigned eqnum = 0;
};


// symbolic equation
class EquationSym {
public:
  unsigned operator[](unsigned i) const {
    return _eq[i].first;
  }

  unsigned &operator[](unsigned i) {
    return _eq[i].first;
  }

  [[nodiscard]] GiNaC::ex coeff(unsigned i) const {
    return _eq[i].second;
  }

  // insert a new item
  void insert(unsigned integral, GiNaC::ex& coeff) {
    _eq.emplace_back(integral, coeff);
  }

  // number of items
  [[nodiscard]] unsigned size() const {
    return _eq.size();
  }

  // get the first integral number
  [[nodiscard]] unsigned first_integral() const {
    return _eq[0].first;
  }

  // get the first coefficient
  GiNaC::ex first_coeff() {
    return _eq[0].second;
  }

  // is empty
  [[nodiscard]] bool empty() const {
    return _eq.empty();
  }

  void sort() {
    std::sort(_eq.begin(), _eq.end(),
              [](const std::pair<unsigned, GiNaC::ex>& i, const std::pair<unsigned, GiNaC::ex>& j){
                return i.first > j.first;
              });
  }

  // clear zero items
  void erase_zero() {
    std::erase_if(_eq, [](const auto &item) { return item.second == 0; });
  }

  // normalize: set the first coeff to one
  void normalize() {
    GiNaC::ex scale = _eq[0].second;
    for (auto &item: _eq)
      item.second /= scale;
  }

  // gauss elimination: eliminate the other equation from this one
  // it is assumed that the other equation has been normalized
  void eliminate(const EquationSym &other, unsigned index) {
    std::vector<std::pair<unsigned, GiNaC::ex>> eq;

    GiNaC::ex scale = _eq[index].second;
    unsigned iother = 0, ithis = 0;
    while (iother < other.size() && ithis < _eq.size()) {
      if(other[iother] > _eq[ithis].first) {
        eq.emplace_back(other[iother], -scale * other._eq[iother].second);
        ++iother;
      }
      else if (other[iother] == _eq[ithis].first) {
        GiNaC::ex coeff = (_eq[ithis].second - scale * other.coeff(iother)).expand();
        if (coeff != 0)
          eq.emplace_back(_eq[ithis].first, coeff);
        ++iother;
        ++ithis;
      }
      else {
        eq.emplace_back(_eq[ithis].first, _eq[ithis].second);
        ++ithis;
      }
    }

    if (iother < other.size())
      for(; iother < other.size(); ++iother)
        eq.emplace_back(other[iother], -scale * other.coeff(iother));
    if (ithis < _eq.size())
      for (; ithis < _eq.size(); ++ithis)
        eq.emplace_back(_eq[ithis]);

    std::swap(eq, _eq);
  }

  // get the underline eq, only readable
  const auto &eq() {
    return _eq;
  }

  // order of equations
  bool operator<(const EquationSym& other) const {
    if (this->first_integral() != other.first_integral())
      return this->first_integral() < other.first_integral();
    else {
      if (this->size() != other.size())
        return this->size() < other.size();
      else
        return this->eqnum < other.eqnum;
    }
  }

  friend std::ostream& operator<<(std::ostream& out, const EquationSym& eq) {
    out << eq._eq;
    return out;
  }

public:
  std::vector<std::pair<unsigned, GiNaC::ex>> _eq;
public:
  unsigned eqnum = 0;
};