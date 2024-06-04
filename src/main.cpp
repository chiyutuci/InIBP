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
  MPRational a1(1, 2);
  MPRational a2(1, 3);
  MPRational a3(1, 4);

  MPRational b1(1, 5);
  MPRational b2(1, 6);
  MPRational b3(1, 7);

  Session session;
  unsigned id = session.new_graph();
  Graph *graph = session.graph(id);

  typedef EvalRationalNumbersData Data;

  std::unique_ptr<EvalRationalNumbers> algptr1(new EvalRationalNumbers());
  std::unique_ptr<Data> data1(new Data());
  EvalRationalNumbers &alg1 = *algptr1;

  std::unique_ptr<EvalRationalNumbers> algptr2(new EvalRationalNumbers());
  std::unique_ptr<Data> data2(new Data());
  EvalRationalNumbers &alg2 = *algptr1;

  std::vector<MPRational> eq1 = {a1, a2, a3};
  std::vector<MPRational> eq2 = {b1, b2, b3};

  alg1.init(std::move(eq1), *data1);
  alg2.init(std::move(eq2), *data2);

  unsigned eq1_id = graph->new_node(std::move(algptr1), std::move(data1), {});
  unsigned eq2_id = graph->new_node(std::move(algptr2), std::move(data2), {});
  std::vector<unsigned> input = {eq1_id, eq2_id};

  std::unique_ptr<Add> addptr(new Add());
  Add &add = *addptr;
  add.init(2, 2);

  unsigned add_id = graph->new_node(std::move(addptr), nullptr, input.data());
  graph->set_output_node(add_id);

  ReconstructionOptions opt;
  opt.max_primes = 20;

  std::unique_ptr<MPRational[]> res(new MPRational[2]);
  Ret ret = session.reconstruct_numeric(id, res.get(), opt);

  std::cout << res[0] << " " << res[1] << std::endl;
  // CLI::App app;

  // std::string configPath;
  // app.add_option("CONFIG", configPath, "The config file to read")
  //     ->default_str("/home/chiyutuci/Works/inibp/example/boxS.yaml")
  //     ->type_name("");

  // CLI11_PARSE(app, argc, argv)

  // std::cout
  //     << "\n--------------- \033[35mInIBP is not IBP\033[0m ---------------"
  //     << std::endl;
  // try {
  //   YAML::Node config =
  //       YAML::LoadFile("/home/chiyutuci/Works/inibp/example/boxS.yaml");
  //   InIBP inibp(config);

  //   inibp.run();
  // } catch (std::runtime_error &e) {
  //   std::cerr << "\n \033[1m\033[31mError:\033[0m " << e.what() << "\n"
  //             << std::endl;
  //   return EXIT_FAILURE;
  // }

  // std::cout << std::endl;
  return EXIT_SUCCESS;
}