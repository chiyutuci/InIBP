#include "family.h"

GiNaC::symtab Family::symtab;

Family::Family(const YAML::Node &config) {
  if (!config["family"])
    throw std::runtime_error("family not found");
  YAML::Node familyConfig = config["family"];

  // check config file keys
  if (!familyConfig["name"])
    throw std::runtime_error("family name not found");
  if (!familyConfig["internals"])
    throw std::runtime_error("family loop momenta not found");
  if (!familyConfig["externals"])
    throw std::runtime_error("family external momenta not found");
  if (!familyConfig["propagators"])
    throw std::runtime_error("family propagators not found");
  if (!familyConfig["sps_rules"])
    throw std::runtime_error("family scalar product rules not found");

  // family name
  _name = familyConfig["name"].as<std::string>();
  // dimension of the family
  if (familyConfig["dimension"]) {
    GiNaC::parser parser;
    _dimension = parser(familyConfig["dimension"].as<std::string>());
  }
  else
    _dimension = GiNaC::possymbol("D");
  if (!is_a<GiNaC::numeric>(_dimension)) {
    _dimension = GiNaC::possymbol("D");
    _symbols.emplace_back("D");
    symtab["D"] = _dimension;
  }
  // loop momenta
  for (const auto &sym: familyConfig["internals"].as<std::vector<std::string>>()) {
    if (symtab.contains(sym))
      throw std::runtime_error("symbol " + sym + " defined more than once");
    _internals.emplace_back(sym);
    symtab[sym] = _internals.back();
  }
  // external momenta
  for (const auto &sym: familyConfig["externals"].as<std::vector<std::string>>()) {
    if (symtab.contains(sym))
      throw std::runtime_error("symbol " + sym + " defined more than once");
    _externals.emplace_back(sym);
    symtab[sym] = _externals.back();
  }
  // invariants
  if (familyConfig["invariants"] && !familyConfig["invariants"].IsNull()) {
    std::string one;
    if (familyConfig["invar_one"] && !familyConfig["invar_one"].IsNull())
      one = familyConfig["invar_one"].as<std::string>();
    for (const auto &inv: familyConfig["invariants"].as<std::vector<std::pair<std::string, unsigned>>>()) {
      if (symtab.contains(inv.first))
        throw std::runtime_error("symbol " + inv.first + " defined more than once");
      _invariants.emplace_back(inv.first, inv.second);
      symtab[inv.first] = _invariants.back().first;
      if (inv.first != one)
        _symbols.emplace_back(inv.first);
    }
  }
  // the invariant set to one
  if (familyConfig["invar_one"] && !familyConfig["invar_one"].IsNull()) {
    auto one = familyConfig["invar_one"].as<std::string>();
    if (!symtab.contains(one))
      throw std::runtime_error("the invariant set to one is not valid");
    _one.append(symtab[one] == 1);
  }

  GiNaC::parser reader(symtab);

  _nints = _internals.size();
  _nexts = _externals.size();
  _nprops = _nints * _nexts + _nints * (_nints + 1) / 2;
  unsigned nsps = _nexts * (_nexts + 1) / 2;

  // invariant scalar products
  if (familyConfig["sps_rules"].size() != nsps)
    throw std::runtime_error("number of scalar product rules does not match the topology");
  for (const auto &sp: familyConfig["sps_rules"].as<std::vector<std::vector<std::string>>>()) {
    if (sp.size() != 3 || !symtab.contains(sp[0]) || !symtab.contains(sp[1]))
      throw std::runtime_error("scalar product rules contains invalid expressions");
    _spsRules.append(reader(sp[0]) * reader(sp[1]) == reader(sp[2]).subs(_one));
  }

  // propagators
  if (familyConfig["propagators"].size() != _nprops)
    throw std::runtime_error("number of propagators does not match the topology");
  for (const auto &prop: familyConfig["propagators"].as<std::vector<std::pair<std::string, std::string>>>())
    _propagators.emplace_back(
            (pow(reader(prop.first), 2) - pow(reader(prop.second), 2)).expand()
                    .subs(_spsRules, GiNaC::subs_options::algebraic).subs(_one).expand());

  std::cout << "\n \033[1m\033[32m#0.0\033[0m   Parsing config file finished." << std::endl;

  _symIndices = generate_symbols("a", _nprops);
  _symProps = generate_symbols("D", _nprops);
}

void Family::init() {
  std::cout << "\n \033[33m#0.1\033[0m   Initializing integral family...\n";
  _compute_sps();
  _compute_symanzik();
  std::cout << "\n \033[1m\033[32m#0.1\033[0m   Initializing integral family finished.\n";

  std::cout << "\n \033[33m#0.2\033[0m   Generating IBP relations...\n";
  _generate_ibp();
  std::cout << "\n \033[1m\033[32m#0.2\033[0m   Generating IBP relations finished.\n";
}

void Family::init_reduce(const YAML::Node &config, Reduce &reduce) const {
  std::cout << "\n \033[33m#0.3\033[0m   Collecting target integrals...\n";
  if (!config["targets"])
    throw std::runtime_error("reduce targets not found");
  reduce._rawTargets = config["targets"].as<std::vector<std::vector<int>>>();

  // find the top sector
  reduce._nprops = _nprops;
  std::vector<bool> hasLine(_nprops, false);
  for (const auto &integral: reduce._rawTargets) {
    if (integral.size() != _nprops) {
      std::stringstream ss;
      ss << integral;
      throw std::runtime_error("target integral " + ss.str() + " is invalid");
    }
    for (unsigned i = 0; i < _nprops; ++i) {
      if (integral[i] > 0 && !hasLine[i])
        hasLine[i] = true;
    }
  }
  for (unsigned i = 0; i < _nprops; ++i)
    if (hasLine[i])
      reduce._top |= 1 << i;
  std::cout << "\n \033[1m\033[32m#0.3\033[0m   Collecting target integrals finished.\n";

  std::cout << "\n \033[33m#0.4\033[0m   Searching trivial sectors...\n";
  _search_trivial_sectors(reduce);
  std::cout << "\n \033[1m\033[32m#0.4\033[0m   Searching trivial sectors finished.\n";

  reduce.prepare_sectors();
}

void Family::print() const {
  std::cout << "\n----------------- \033[36mFamily Info\033[0m ------------------\n"
            << "\n  Topology: " << _name << "   Dimension: " << _dimension
            << "\n\n  Internals:";
  for (const auto &sym: _internals)
    std::cout << " " << sym;

  std::cout << "\n\n  Externals:";
  for (const auto &sym: _externals)
    std::cout << " " << sym;

  std::cout << "\n\n  Invariants:";
  for (const auto &inv: _invariants)
    std::cout << " " << inv.first;

  std::cout << "\n\n  Propagators: " << _nprops << "   IBP relations: " << _ibp.size();

  std::cout << std::endl;
}

void Family::_generate_ibp() {
  GiNaC::ex coeff, coeffD;
  std::map<std::vector<int>, GiNaC::ex> equation;

  // i: l_i in derivatives
  // j: l_j or p_j in nominators
  // s: prop_s in denominators
  for (size_t i = 0; i < _nints; ++i) {
    for (size_t j = 0; j < _nints + _nexts; ++j) {
      equation.clear();
      // cases of i == j generate g^u_u = D
      if (i == j)
        equation[std::vector<int>(_nprops, 0)] += _dimension;
      for (size_t s = 0; s < _nprops; ++s) {
        if (j < _nints)
          coeff = (-_symIndices[s]) * _internals[j] * GiNaC::diff(_propagators[s], _internals[i]);
        else
          coeff = (-_symIndices[s]) * _externals[j - _nints] *
                  GiNaC::diff(_propagators[s], _internals[i]);
        coeff = coeff.expand().subs(_spsRules, GiNaC::subs_options::algebraic).expand();
        if (coeff != 0) {
          std::vector<int> integral(_nprops, 0);
          integral[s] = 1;
          // substitute scalar products by propagators
          coeff = coeff.subs(_spsFromProps, GiNaC::subs_options::algebraic);
          // t: D_t in coeff
          for (size_t t = 0; t < _nprops; ++t) {
            coeffD = coeff.diff(_symProps[t]);
            if (coeffD != 0) {
              integral[t] -= 1;
              equation[integral] += coeffD;
              integral[t] += 1;
            }
            coeff = coeff.subs(_symProps[t] == 0);
            if (coeff == 0)
              break;
          }
          coeff = coeff.expand();
          if (coeff != 0)
            equation[integral] += coeff;
        }
      }
      // add to _ibp
      IBPProto ibp;
      for (const auto &item: equation | std::ranges::views::reverse) {
        coeff = item.second.expand();
        if (coeff != 0)
          ibp.emplace_back(item.first, coeff);
      }
      if (!ibp.empty())
        _ibp.push_back(std::move(ibp));
    }
  }
}

void Family::_search_trivial_sectors(Reduce &reduce) const {
  GiNaC::ex gPoly = _uPoly + _fPoly;
  GiNaC::ex gDiff;

  std::vector<GiNaC::possymbol> kvec = generate_symbols("k", _nprops);
  GiNaC::lst klst;
  for (const auto &sym: kvec)
    klst.append(sym);

  for (unsigned i = 0; i < _nprops; ++i)
    gDiff += kvec[i] * _symIndices[i] * gPoly.diff(_symIndices[i]);
  gDiff = (gDiff - gPoly).expand();

  reduce._sectors = std::vector<bool>(reduce._top + 1, false);
  std::vector<unsigned> trivials;
  for (unsigned sector = reduce._top; sector >= (1 << _nints) - 1; --sector) {
    // check if this sector belongs to top sector
    if ((sector & reduce._top) == sector) {
      // check if this sector belongs to a trivial sector
      if (std::find_if(trivials.begin(), trivials.end(), [sector](unsigned trivial) {
        return (sector & trivial) == sector;
      }) != trivials.end())
        continue;

      GiNaC::lst zeros;
      for (unsigned i = 0; i < _nprops; ++i)
        if ((sector & (1 << i)) == 0)
          zeros.append(_symIndices[i] == 0);
      GiNaC::ex gSector = gDiff.subs(zeros);

      // collect the k equations
      std::map<std::vector<unsigned>, GiNaC::ex> coeffs;
      for (unsigned item = 0; item < gSector.nops(); ++item) {
        GiNaC::ex term = gSector.op(item);
        std::vector<unsigned> indices(_nprops, 0);
        for (unsigned i = 0; i < _nprops; ++i) {
          if (sector & (1 << i)) {
            if (term.coeff(_symIndices[i], 2) != 0) {
              term = term.coeff(_symIndices[i], 2);
              indices[i] = 2;
            }
            else if (term.coeff(_symIndices[i], 1) != 0) {
              term = term.coeff(_symIndices[i], 1);
              indices[i] = 1;
            }
          }
        }
        coeffs[indices] += term;
      }
      GiNaC::lst eqs;
      for (const auto &item: coeffs) {
        GiNaC::ex coeff = item.second.expand();
        if (coeff != 0)
          eqs.append(coeff == 0);
      }
      // solve the k equations
      if (nops(GiNaC::lsolve(eqs, klst)) == 0)
        reduce._sectors[sector] = true;
      else
        trivials.push_back(sector);
    }
  }
}

void Family::_compute_sps() {
  GiNaC::matrix propSpMat(_nprops, _nprops);
  GiNaC::matrix propSpConst(_nprops, 1);

  // set all variables to zero to get constant items
  GiNaC::lst allZero;
  for (unsigned i = 0; i < _nints; ++i) {
    for (unsigned j = i; j < _nints; ++j)
      allZero.append(_internals[i] * _internals[j] == 0);
    for (unsigned j = 0; j < _nexts; ++j)
      allZero.append(_internals[i] * _externals[j] == 0);
  }
  // compute the matrix from sp to prop
  for (unsigned s = 0; s < _nprops; ++s) {
    propSpConst(s, 0) = _propagators[s].subs(allZero, GiNaC::subs_options::algebraic);

    unsigned index = 0;
    for (unsigned i = 0; i < _nints; ++i) {
      GiNaC::ex coeff = _propagators[s].diff(_internals[i]);
      propSpMat(s, index) = coeff.diff(_internals[i]) / 2;
      ++index;
      for (unsigned j = i + 1; j < _nints; ++j) {
        propSpMat(s, index) = coeff.diff(_internals[j]);
        ++index;
      }
      for (unsigned j = 0; j < _nexts; ++j) {
        propSpMat(s, index) = coeff.diff(_externals[j]);
        ++index;
      }
    }
  }
  // check if the matrix is invertible
  if (propSpMat.rank() < _nprops)
    throw std::runtime_error("the set of propagators is incomplete");
  // compute the expression for sps
  GiNaC::matrix spPropVec(_nprops, 1);
  for (unsigned i = 0; i < _nprops; ++i)
    spPropVec(i, 0) = _symProps[i] - propSpConst(i, 0);
  GiNaC::matrix propSpInv = propSpMat.inverse();
  spPropVec = propSpInv.mul(spPropVec);
  // fill the rule list
  unsigned index = 0;
  for (unsigned i = 0; i < _nints; ++i) {
    for (unsigned j = i; j < _nints; ++j) {
      _spsFromProps.append(_internals[i] * _internals[j] == spPropVec(index, 0));
      ++index;
    }
    for (unsigned j = 0; j < _nexts; ++j) {
      _spsFromProps.append(_internals[i] * _externals[j] == spPropVec(index, 0));
      ++index;
    }
  }
}

void Family::_compute_symanzik() {
  GiNaC::ex schwinger;
  for (unsigned i = 0; i < _nprops; ++i)
    schwinger += -_symIndices[i] * _propagators[i];

  // the constant item
  GiNaC::lst allZero;
  for (const auto &ex: _internals)
    allZero.append(ex == 0);
  GiNaC::ex J = -schwinger.subs(allZero, GiNaC::subs_options::algebraic);

  // matrix l.M.l
  GiNaC::matrix M(_nints, _nints);
  // vector l.v
  GiNaC::matrix V(_nints, 1);

  for (unsigned i = 0; i < _nints; ++i) {
    GiNaC::ex di = schwinger.diff(_internals[i]);
    for (unsigned j = i; j < _nints; ++j) {
      M(i, j) = di.diff(_internals[j]) / 2;
      if (i != j)
        M(j, i) = M(i, j);
    }
    V(i, 0) = -di.subs(allZero, GiNaC::subs_options::algebraic) / 2;
  }

  _uPoly = M.determinant().expand();
  _fPoly = (_uPoly * J + V.transpose().mul(M.inverse().mul_scalar(_uPoly)).mul(V)(0, 0)).expand()
          .subs(_spsRules, GiNaC::subs_options::algebraic).expand();
}

std::vector<GiNaC::possymbol> Family::generate_symbols(const std::string &name, unsigned n) {
  std::vector<GiNaC::possymbol> symbols(n);
  for (size_t i = 0; i < n; ++i)
    symbols[i] = GiNaC::possymbol(name + std::to_string(i + 1));
  return symbols;
}

void Reduce::print() const {
  std::cout << "\n----------------- \033[36mReduce Info\033[0m ------------------\n"
            << "\n  Targets: " << _rawTargets.size() << "   Top sector: " << _top;

  std::cout << "\n\n  Non-trivial sectors: "
            << std::count_if(_sectors.begin(), _sectors.end(), [](bool value) {
              return value;
            });

  std::cout << std::endl;
}

void Reduce::prepare_sectors() {
  _lines = std::vector<bool>(_nprops, false);
  for (unsigned i = 0; i < _nprops; ++i)
    if (_top & (1 << i))
      _lines[i] = true;

  // sort non-trivial sectors
  unsigned nsec = std::count_if(_sectors.begin(), _sectors.end(), [](bool value) {
    return value;
  });
  std::vector<unsigned> sectors(nsec);
  for (unsigned i = 0, j = nsec - 1; i < _sectors.size(); ++i)
    if (_sectors[i])
      sectors[j--] = i;
  std::sort(sectors.begin(), sectors.end(), [](unsigned a, unsigned b) {
    unsigned alines = std::popcount(a);
    unsigned blines = std::popcount(b);
    if (alines > blines)
      return true;
    else if (alines == blines)
      return a > b;
    else
      return false;
  });

  // fill the reduce sectors
  _reduceSectors = std::vector<Sector>(nsec);
  for (unsigned i = 0; i < sectors.size(); ++i) {
    _reduceSectors[i]._id = sectors[i];
    for (unsigned j = 0; j < _nprops; ++j)
      if ((sectors[i] & (1 << j)) == 0 && _lines[j])
        _reduceSectors[i]._superSectors.push_back(sectors[i] | (1 << j));
    for (unsigned j = 0; j < _nprops; ++j)
      if (sectors[i] & (1 << j) && _sectors[sectors[i] ^ (1 << j)])
        _reduceSectors[i]._subSectors.push_back(sectors[i] ^ (1 << j));
  }
}
