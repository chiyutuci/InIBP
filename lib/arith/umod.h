#pragma once

#include <iostream>

#include "flint/nmod_vec.h"


typedef unsigned long uint64;
typedef signed long sint64;

const nmod_t MOD64{9223372036854775783, 50, 1};

// finite field over 9223372036854775783 (the largest 63 bit prime)
// use umod64::from() to convert a signed integer to umod64
class umod64 {
public:
  umod64() = default;

  // it is assumed that num is in the range [0, MOD64)
  explicit umod64(uint64 num) : _num(num) {}

  friend umod64 operator-(const umod64 &num);
  friend std::ostream &operator<<(std::ostream &out, const umod64 &num);

  // convert a signed integer to a umod64
  static umod64 from(sint64 num) {
    umod64 res;
    if (num > 0) {
      NMOD_RED(res._num, num, MOD64);
    }
    else if (num < 0) {
      NMOD_RED(res._num, -num, MOD64);
      return -res;
    }
    return res;
  }

  // operators

  umod64 operator+(const umod64 &other) const {
    return umod64{_nmod_add(_num, other._num, MOD64)};
  }

  umod64 operator-(const umod64 &other) const {
    return umod64{_nmod_sub(_num, other._num, MOD64)};
  }

  umod64 operator*(const umod64 &other) const {
    return umod64{nmod_mul(_num, other._num, MOD64)};
  }

  umod64 operator/(const umod64 &other) const {
    return umod64{nmod_div(_num, other._num, MOD64)};
  }

  umod64 operator^(uint64 pow) const {
    return umod64{nmod_pow_ui(_num, pow, MOD64)};
  }

  umod64 &operator+=(const umod64 &other) {
    _num = _nmod_add(_num, other._num, MOD64);
    return *this;
  }

  umod64 &operator-=(const umod64 &other) {
    _num = _nmod_sub(_num, other._num, MOD64);
    return *this;
  }

  umod64 &operator*=(const umod64 &other) {
    _num = nmod_mul(_num, other._num, MOD64);
    return *this;
  }

  umod64 &operator/=(const umod64 &other) {
    _num = nmod_div(_num, other._num, MOD64);
    return *this;
  }

  bool operator==(const umod64 &other) const { return _num == other._num; }

  bool operator!=(const umod64 &other) const { return _num != other._num; }

private:
  uint64 _num = 0;
};

// negative
inline umod64 operator-(const umod64 &num) {
  return umod64{nmod_neg(num._num, MOD64)};
}

std::ostream &operator<<(std::ostream &out, const umod64 &num) {
  out << num._num;
  return out;
}