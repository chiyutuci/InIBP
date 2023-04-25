#include "inibp.h"

InIBP::InIBP(const YAML::Node &node) : _family(node) {
	try {
		_reduce = node["reduce"].as<Reduce>();
	}
	catch (YAML::Exception &e) {
		std::cerr << "InIBP: " << e.what() << std::endl;
	}
}

void InIBP::init() {

	_family.print();
	_reduce.print();
}

void InIBP::run() {}


