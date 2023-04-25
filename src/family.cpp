#include "family.h"

GiNaC::symtab Family::symtab;

Family::Family(const YAML::Node &config) {
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
	if (!familyConfig["sp_rules"])
		throw std::runtime_error("family scalar product rules not found");

	// family name
	_name = familyConfig["name"].as<std::string>();
	// dimension of the family
	if (familyConfig["dimension"]) {
		GiNaC::parser parser;
		_dimension = parser(familyConfig["dimension"].as<std::string>());
	} else
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
	if (familyConfig["invariants"]) {
		std::string one;
		if (familyConfig["invar_one"])
			one = familyConfig["invar_one"].as<std::string>();
		for (const auto &inv: familyConfig["invariants"]
						.as<std::vector<std::pair<std::string, unsigned>>>()) {
			if (symtab.contains(inv.first))
				throw std::runtime_error("symbol " + inv.first +
				                         " defined more than once");
			_invariants.emplace_back(inv.first, inv.second);
			symtab[inv.first] = _invariants.back().first;
			if (inv.first != one)
				_symbols.emplace_back(inv.first);
		}
	}
	// the invariant set to one
	if (familyConfig["invar_one"]) {
		auto one = familyConfig["invar_one"].as<std::string>();
		if (!symtab.contains(one))
			throw std::runtime_error("the invariant set to one is not valid");
		_one = (symtab[one] == 1);
	}

	GiNaC::parser reader(symtab);

	_nints = _internals.size();
	_nexts = _externals.size();
	_nprops = _nints * (_nexts + (_nints + 1) / 2);
	unsigned nsps = _nexts * (_nexts + 1) / 2;

	// invariant scalar products
	if (familyConfig["sp_rules"]) {
		if (familyConfig["sp_rules"].size() != nsps)
			throw std::runtime_error("number of scalar product rules does not match the topology");
		for (const auto &sp: familyConfig["sp_rules"]
						.as<std::vector<std::vector<std::string>>>()) {
			if (!symtab.contains(sp[0]) || !symtab.contains(sp[1]))
				throw std::runtime_error("scalar product rules contains invalid expressions");
			_spsRules.append(reader(sp[0]) * reader(sp[1]) == reader(sp[2]).subs(_one).expand());
		}
	}
	// propagators
	if (familyConfig["propagators"].size() != _nprops)
		throw std::runtime_error("number of propagators does not match the topology");
	for (const auto &prop: familyConfig["propagators"]
					.as<std::vector<std::pair<std::string, std::string>>>())
		_propagators.emplace_back(
						(pow(reader(prop.first), 2) - pow(reader(prop.second), 2)).subs(_one).expand()
										.subs(_spsRules, GiNaC::subs_options::algebraic)).expand();

	_symIndices = generate_symbols("a", _nprops);
	_symProps = generate_symbols("D", _nprops);

	_compute_sps();
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

void Family::generate_ibp() {
	std::cout << "\n \033[32m#0.1\033[0m   Generating IBP relations..." << std::endl;

	_generate_ibp_detail();

	std::cout << "\n \033[32m#0.1\033[0m   Generating IBP relations finished." << std::endl;
}

void Family::_generate_ibp_detail() {
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
					coeff = (-_symIndices[s]) * _internals[j]
					        * GiNaC::diff(_propagators[s], _internals[i]);
				else
					coeff = (-_symIndices[s]) * _externals[j - _nints] *
					        GiNaC::diff(_propagators[s], _internals[i]);
				coeff = coeff.expand().subs(_spsRules, GiNaC::subs_options::algebraic).expand();
				if (coeff != 0) {
					std::vector<int> integral(_nprops, 0);
					integral[s] = 1;
					// substitute scalar products by propagators
					coeff = coeff.expand().subs(_spsFromProps, GiNaC::subs_options::algebraic);
					// t: D_t in coeff
					for (size_t t = 0; t < _nprops; ++t) {
						coeffD = coeff.diff(_symProps[t]);
						if (coeffD != 0) {
							integral[t] -= 1;
							equation[integral] += coeffD;
							integral[t] += 1;
						}
						coeff = coeff.subs(_symProps[t] == 0);
						if (coeff == 0) break;
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
			_spsFromProps.append(_internals[i] * _externals[j] ==
			                     spPropVec(index, 0));
			++index;
		}
	}
}

std::vector<GiNaC::possymbol> Family::generate_symbols(const std::string &name, unsigned n) {
	std::vector<GiNaC::possymbol> symbols(n);
	for (size_t i = 0; i < n; ++i)
		symbols[i] = GiNaC::possymbol(name + std::to_string(i + 1));
	return symbols;
}

Reduce::Reduce(const YAML::Node &node) {
	YAML::Node reduceConfig = node["reduce"];

	_top = reduceConfig["top"].as<unsigned>();
	_posi = reduceConfig["posi"].as<unsigned>();
	_rank = reduceConfig["rank"].as<unsigned>();
	_dot = reduceConfig["dot"].as<unsigned>();
}

void Reduce::print() const {
	std::cout << "\n----------------- \033[36mReduce Info\033[0m ------------------\n"
	          << "\n  Top: " << _top << "   Posi: " << _posi << "   Rank: " << _rank
	          << "   Dot: " << _dot << std::endl;
}
