#pragma once

#include <string>
#include <ranges>
#include <iostream>

#include "ginac/ginac.h"
#include "yaml-cpp/yaml.h"

#include "utils.h"
#include "sector.h"


class Reduce;


class Family {
public:
  explicit Family(const YAML::Node &config);

  // check if the reduction job is valid
  void check_reduce(const Reduce &) const;

  // initialize the family
  void init();
  // prepare reduce information
  void init_reduce(Reduce &) const;
  // print family information
  void print() const;

  static std::vector<GiNaC::possymbol> generate_symbols(const std::string &, unsigned);

private:
  // compute expressions of scalar products over propagators
  void _compute_sps();
  //compute symanzik polynomials
  void _compute_symanzik();

  // generate ibp relations
  void _generate_ibp();
  // search trivial sectors
  void _search_trivial_sectors(Reduce &) const;

public:
  static GiNaC::symtab symtab;

private:
  // family name
  std::string _name;
  // dimension of the family
  GiNaC::ex _dimension;
  // loop momenta
  std::vector<GiNaC::possymbol> _internals;
  // external momenta
  std::vector<GiNaC::possymbol> _externals;
  // number of internals
  unsigned _nints = 0;
  // number of externals
  unsigned _nexts = 0;
  // kinematics invariants
  //  first: name
  //  second: dimension
  std::vector<std::pair<GiNaC::possymbol, unsigned>> _invariants;
  // symbols
  std::vector<GiNaC::possymbol> _symbols;
  // propagators
  std::vector<GiNaC::ex> _propagators;
  // scalar products rules
  GiNaC::lst _spsRules;
  // scalar products over propagators
  GiNaC::lst _spsFromProps;
  // the invariant set to one
  GiNaC::lst _one;
  // number of propagators
  unsigned _nprops = 0;

  // symbols for integral indices
  std::vector<GiNaC::possymbol> _symIndices;
  // symbols for integral propagators
  std::vector<GiNaC::possymbol> _symProps;
  // symanzik U
  GiNaC::ex _uPoly;
  // symanzik F
  GiNaC::ex _fPoly;
  // ibp relation
  // first:  integral indices
  // second: coefficient
  typedef std::vector<std::pair<std::vector<int>, GiNaC::ex>> IBPProto;
  // ibp relations prototype
  std::vector<IBPProto> _ibp;

};


class Reduce {
public:
  explicit Reduce(const YAML::Node &config);

  // print reduciton job info
  void print() const;


  friend class Family;


private:
  // top level sector
  unsigned _top = 0;
  // sum of positive indices
  unsigned _posi = 0;
  // sum of negative indices
  unsigned _rank = 0;
  // posi - lines
  unsigned _dot = 0;
  // non-trivial sectors
  // true: non-trivial
  // false: trivial or irrelevant
  std::vector<bool> _sectors;
  // the reduction jobs
  std::vector<Sector> _reduceSectors;
};