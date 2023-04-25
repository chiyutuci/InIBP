#pragma once

#include <fstream>
#include <iostream>

#include <yaml-cpp/yaml.h>

#include "family.h"

class InIBP {
public:
	explicit InIBP(const YAML::Node &);

public:
	// prepare family and reduction job
	void init();

	// run the reduction
	void run();

private:
private:
	Family _family;
	Reduce _reduce;
};