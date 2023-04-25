#include "inibp.h"

InIBP::InIBP(const YAML::Node &node) : _family(node), _reduce(node) {
	std::cout << "\n \033[32m#0.0\033[0m   Parsing config file finished." << std::endl;
}

void InIBP::init() {

	_family.print();
	_reduce.print();
}

void InIBP::run() {}


