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
	} else _dimension = GiNaC::possymbol("D");
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
		for (const auto &inv: familyConfig["invariants"].
						as<std::vector<std::pair<std::string, unsigned >>>()) {
			if (symtab.contains(inv.first))
				throw std::runtime_error("symbol " + inv.first + " defined more than once");
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

	unsigned nint = _internals.size();
	unsigned next = _externals.size();
	unsigned nsps = next * (next + 1) / 2;
	unsigned nprops = nint * (next + (nint + 1) / 2);
	// invariant scalar products
	if (familyConfig["sp_rules"]) {
		if (familyConfig["sp_rules"].size() != nsps)
			throw std::runtime_error("number of scalar product rules does not match the topology");
		for (const auto &sp: familyConfig["sp_rules"].
						as<std::vector<std::vector<std::string>>>()) {
			if (!symtab.contains(sp[0]) || !symtab.contains(sp[1]))
				throw std::runtime_error("scalar product rules contains invalid expressions");
			_sp_rules.append(reader(sp[0]) * reader(sp[1]) == reader(sp[2]).subs(_one).expand());
		}
	}
	// propagators
	if (familyConfig["propagators"].size() != nprops)
		throw std::runtime_error("number of propagators does not match the topology");
	for (const auto &prop: familyConfig["propagators"].
					as<std::vector<std::pair<std::string, std::string>>>())
		_propagators.emplace_back(
						(pow(reader(prop.first), 2) - pow(reader(prop.second), 2)).subs(_one).expand()
										.subs(_sp_rules, GiNaC::subs_options::algebraic)).expand();
}

void Family::print() const {
	std::cout << "\n----------------- \033[36mFamily Info\033[0m ------------------"
	          << "\n  Topology: " << _name << "   Dimension: " << _dimension
	          << "\n  Internals:";
	for (const auto &sym: _internals)
		std::cout << " " << sym;

	std::cout << "\n  Externals:";
	for (const auto &sym: _externals)
		std::cout << " " << sym;

	std::cout << "\n  Invariants:";
	for (const auto &inv: _invariants)
		std::cout << " " << inv.first;

	std::cout << std::endl;
}


Reduce::Reduce(const YAML::Node &node) {
	YAML::Node reduceConfig = node["reduce"];

	_top = reduceConfig["top"].as<unsigned>();
	_posi = reduceConfig["posi"].as<unsigned>();
	_rank = reduceConfig["rank"].as<unsigned>();
	_dot = reduceConfig["dot"].as<unsigned>();
}

void Reduce::print() const {
	std::cout << "\n----------------- \033[36mReduce Info\033[0m ------------------"
	          << "\n  Top: " << _top << "   Posi: " << _posi
	          << "   Rank: " << _rank << "   Dot: " << _dot
	          << std::endl;
}
