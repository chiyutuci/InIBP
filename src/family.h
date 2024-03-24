#pragma once

#include <string>
#include <ranges>
#include <random>
#include <iostream>

#include "ginac/ginac.h"
#include "yaml-cpp/yaml.h"

#include "utils.h"
#include "sector.h"


class Family {
public:
  explicit Family(const YAML::Node &);

  // initialize the family
  void init();
  // prepare reduce information
  void init_reduce(const YAML::Node &, Reduce &) const;
  // print family information
  void print() const;

  // run the reduction
  void run_reduce(Reduce &) const;

  static std::vector<GiNaC::symbol> generate_symbols(const std::string &, unsigned);

private:
  // compute expressions of scalar products over propagators
  void _compute_sps();
  //compute symanzik polynomials
  void _compute_symanzik();

  // generate ibp relations
  void _generate_ibp();
  // generate li relations
  void _generate_li();
  // generate ibp over finite filed
  void _generate_ibp_ff();
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
  std::vector<GiNaC::symbol> _symIndices;
  // symbols for integral propagators
  std::vector<GiNaC::symbol> _symProps;
  // symanzik U
  GiNaC::ex _uPoly;
  // symanzik F
  GiNaC::ex _fPoly;
  // ibp relations prototype
  std::vector<IBPProto> _ibp;
  // ibp relations prototype over finite field
  std::vector<IBPProtoFF> _ibpFF;
};

class Reduce {
public:
  // prepare sectors relations
  void prepare_sectors();

  // print reduciton job info
  void print() const;

  friend class Family;

private:
  // symbols
  std::vector<GiNaC::possymbol> _symbols;
  std::vector<GiNaC::symbol> _symIndices;
  // raw target integrals
  std::vector<RawIntegral> _rawTargets;
  // top sector
  unsigned _top = 0;
  // number of propagators
  unsigned _nprops = 0;
  // top sector lines
  std::vector<bool> _lines;
  // non-trivial sectors
  // true: non-trivial
  // false: trivial or irrelevant
  std::vector<bool> _sectors;
  // the reduction jobs
  std::vector<Sector> _reduceSectors;
};