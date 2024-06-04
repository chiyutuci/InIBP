#include <chrono>
#include <ctime>
#include <iostream>

#include "CLI11.hpp"
#include "yaml-cpp/yaml.h"

#include "inibp.h"
#include <flint/nmod_vec.h>

#include <fflow/alg_functions.hh>
#include <fflow/alg_lists.hh>
#include <fflow/algorithm.hh>
#include <fflow/graph.hh>

using namespace fflow;

int main(int argc, char **argv) {
  auto now = std::chrono::system_clock::now();
  std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
  std::tm *localtime = std::localtime(&currentTime);
  std::cout << "\nProgram begin: " << std::asctime(localtime);

  CLI::App app;

  std::string configPath;
  app.add_option("CONFIG", configPath, "The config file to read")
      ->type_name("");

  CLI11_PARSE(app, argc, argv)

  std::cout
      << "\n--------------- \033[35mInIBP is not IBP\033[0m ---------------"
      << std::endl;
  try {
    YAML::Node config = YAML::LoadFile(configPath);
    // YAML::Node config =
    //     YAML::LoadFile("/home/chiyutuci/Works/inibp/example/1.yaml");
    InIBP inibp(config);

    inibp.run();
  } catch (std::runtime_error &e) {
    std::cerr << "\n \033[1m\033[31mError:\033[0m " << e.what() << "\n"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  return EXIT_SUCCESS;
}