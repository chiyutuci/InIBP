#include "CLI11.hpp"
#include <iostream>

#include "inibp.h"

int main(int argc, char **argv) {
  CLI::App app;

  std::string config;
  app.add_option("CONFIG", config, "The config file to read")
      ->type_name("")
      ->required()
      ->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  return 0;
}