#pragma once

#include <iostream>

#include "flint/fmpz.h"
#include "flint/nmod_vec.h"


typedef unsigned long uint64;
typedef signed long sint64;

const nmod_t MOD64{9223372036854775783, 50, 1};
const uint64 PRIMES64[] = {
        8646118801249252579,
        8283289069716051751,
        7708750630245640033,
        7266804057017283853,
        6956220289551338803,
        6379266932664197357,
        5848477568814077963,
        5152842736579902761,
        4650640856291619181,
        4230608376698653057,
        3762356528626792027,
        3219529993405051873,
        2910658703784597323,
        2276848867800693079,
        1925632528975719227,
        1277294106943470761
};

// finite field over 9223372036854775783 (the largest 63 bit prime)
// use umod64::from() to convert a signed integer to umod64
class umod64 {
public:
  umod64() = default;

  // it is assumed that num is in the range [0, MOD64)
  explicit umod64(uint64 num) : _num(num) {}

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

  // convert a string to a umod64
  static umod64 from(const std::string &s) {
    fmpz_t num;
    fmpz_init(num);
    fmpz_set_str(num, s.c_str(), 10);

    umod64 res;
    res._num = fmpz_get_nmod(num, MOD64);

    fmpz_clear(num);
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

  bool operator==(uint64 other) const { return _num == other; }

  bool operator!=(uint64 other) const { return _num != other; }

  friend umod64 operator-(const umod64 &num) {
    return umod64{nmod_neg(num._num, MOD64)};
  }

  friend std::ostream &operator<<(std::ostream &out, const umod64 &num) {
    out << num._num;
    return out;
  }

private:
  uint64 _num = 0;
};