#include <iostream>

#include "CLI11.hpp"
#include "yaml-cpp/yaml.h"

#include "inibp.h"

int main(int argc, char **argv) {
	CLI::App app;

	std::string configPath;
	app.add_option("CONFIG", configPath, "The config file to read")
					->type_name("")
					->required()
					->check(CLI::ExistingFile);

	CLI11_PARSE(app, argc, argv);

	YAML::Node config = YAML::LoadFile(configPath);
	InIBP inibp(config);

	inibp.init();
	inibp.run();

	return 0;
}