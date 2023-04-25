#pragma once

#include <string>

#include "ginac/ginac.h"
#include "yaml-cpp/yaml.h"

class Family {
public:
	explicit Family(const YAML::Node &config);

	// print family info
	void print() const;

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
	// kinematics invariants
	//  first: name
	//  second: dimension
	std::vector<std::pair<GiNaC::possymbol, unsigned >> _invariants;
	// the invariant set to one
	GiNaC::ex _one;
	// symbols
	std::vector<GiNaC::possymbol> _symbols;

	// substitution rules
	GiNaC::lst _sp_rules;
	// propagators
	std::vector<GiNaC::ex> _propagators;
};

class Reduce {
public:
	explicit Reduce(const YAML::Node &config);

	// print reduciton job info
	void print() const;

private:
	// top level sector
	unsigned _top = 0;
	// sum of positive indices
	unsigned _posi = 0;
	// sum of negative indices
	unsigned _rank = 0;
	// posi - lines
	unsigned _dot = 0;
};