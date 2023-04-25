#include "family.h"

GiNaC::symtab Family::symtab;

Family::Family(const YAML::Node &config) {
	try {
		YAML::Node familyConfig = config["family"];

		_name = familyConfig["name"].as<std::string>();

		for (const auto &ex: familyConfig["internals"].as<std::vector<std::string>>()) {
			symtab[ex] = GiNaC::possymbol(ex);
			_internals.emplace_back(symtab[ex]);
		}

		for (const auto &ex: familyConfig["externals"].as<std::vector<std::string>>()) {
			symtab[ex] = GiNaC::possymbol(ex);
			_externals.emplace_back(symtab[ex]);
		}
	}
	catch (YAML::Exception &e) {
		std::cerr << "Family: " << e.what() << std::endl;
	}
}

void Family::print() const {
	std::cout << "\n------------- Family Info -------------"
	          << "\n  Topology: " << _name
	          << "\n  Internal:";
	for (const auto &ex: _internals) {
		std::cout << " " << ex;
	}
	std::cout << "\n  External:";
	for (const auto &ex: _externals) {
		std::cout << " " << ex;
	}
	std::cout << std::endl;
}

void Reduce::print() const {
	std::cout << "\n------------- Reduce Info -------------"
	          << "\n  Top: " << top << "  Posi: " << posi
	          << "  Rank: " << rank << "  Dot: " << dot
	          << std::endl;
}
