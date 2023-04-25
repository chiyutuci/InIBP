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

	std::cout << "\n--------------- \033[34mInIBP is not IBP\033[0m ---------------" << std::endl;
	try {
		YAML::Node config = YAML::LoadFile(configPath);
		InIBP inibp(config);

		inibp.init();
		inibp.run();
	}
	catch (std::runtime_error &e) {
		std::cerr << " \033[31mError:\033[0m " << e.what()
		          << "\n" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << std::endl;
	return EXIT_SUCCESS;
}