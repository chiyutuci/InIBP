#pragma once

#include <string>
#include <ranges>
#include <iostream>

#include "ginac/ginac.h"
#include "yaml-cpp/yaml.h"

#include "utils.h"


class Family {
public:
	explicit Family(const YAML::Node &config);

	// genetate ibp relations
	void generate_ibp();
	// print family info
	void print() const;

	static std::vector<GiNaC::possymbol> generate_symbols(const std::string &, unsigned);

private:
	// generate ibp relations
	void _generate_ibp_detail();
	// compute expressions of scalar products over propagators
	void _compute_sps();

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
	GiNaC::ex _one;
	// number of propagators
	unsigned _nprops = 0;

	// symbols for integral indices
	std::vector<GiNaC::possymbol> _symIndices;
	// symbols for integral propagators
	std::vector<GiNaC::possymbol> _symProps;
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