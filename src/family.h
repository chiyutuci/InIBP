#pragma once

#include <string>

#include "ginac/ginac.h"
#include "yaml-cpp/yaml.h"

class Family {
public:
	explicit Family(const YAML::Node &config);

	void print() const;

public:
	static GiNaC::symtab symtab;

private:
	std::string _name;
	std::vector<GiNaC::ex> _internals;
	std::vector<GiNaC::ex> _externals;
};

struct Reduce {
	// top level sector
	unsigned top = 0;
	// sum of positive indices
	unsigned posi = 0;
	// sum of negative indices
	unsigned rank = 0;
	// posi - lines
	unsigned dot = 0;

	void print() const;
};


namespace YAML {
	template<>
	struct convert<Reduce> {
		static bool decode(const Node &node, Reduce &reduce) {
			if (!node.IsMap())
				return false;

			reduce.top = node["top"].as<unsigned>();
			reduce.posi = node["posi"].as<unsigned>();
			reduce.rank = node["rank"].as<unsigned>();
			reduce.dot = node["dot"].as<unsigned>();

			return true;
		}
	};
} // namespace YAML