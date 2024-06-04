#include "sector.h"

#include <fflow/alg_functions.hh>
#include <fflow/graph.hh>
#include <fflow/numeric_solver.hh>
#include <fstream>

using namespace fflow;

std::map<std::pair<int, int>, std::vector<std::vector<int>>>
    Sector::combinations;

void Sector::generate_combinations(int number, int sum) {
  for (int num = 0; num <= number; ++num)
    for (int s = 0; s <= sum; ++s) {
      if (combinations.contains({num, s}))
        continue;
      else if (num == 0)
        combinations[{num, s}] = {{}};
      else if (num == 1)
        combinations[{num, s}] = {{s}};
      else if (s == 0)
        combinations[{num, s}] = {std::vector<int>(num, 0)};
      else
        for (int i = 0; i <= s; ++i) {
          for (const auto &comb : combinations[{num - 1, s - i}]) {
            std::vector<int> newComb = comb;
            newComb.emplace_back(i);
            combinations[{num, s}].emplace_back(std::move(newComb));
          }
        }
    }
}

void Sector::_generate_seeds() {
  int lines = std::popcount(_id);
  int zeros = (int)_nprops - lines;
  int number = std::max(lines, zeros);
  int sum = std::max((int)_depth - std::popcount(_id), (int)_rank);
  generate_combinations(number, sum);

  for (int depth = 0; depth <= _depth - std::popcount(_id); ++depth) {
    if (zeros != 0)
      for (int rank = 0; rank <= _rank; ++rank)
        for (const auto &depthComb : combinations[{lines, depth}])
          for (const auto &rankComb : combinations[{zeros, rank}]) {
            RawIntegral integral(_lines);
            unsigned posLines = 0;
            unsigned posZeros = 0;
            for (unsigned i = 0; i < _nprops; i++) {
              if (_lines[i])
                integral[i] += depthComb[posLines++];
              else
                integral[i] -= rankComb[posZeros++];
            }
            _seeds.emplace_back(std::move(integral));
          }
    else {
      for (const auto &depthComb : combinations[{lines, depth}]) {
        RawIntegral integral(_lines);
        unsigned posLines = 0;
        for (unsigned i = 0; i < _nprops; i++)
          if (_lines[i])
            integral[i] += depthComb[posLines++];
        _seeds.emplace_back(std::move(integral));
      }
    }
  }
  for (unsigned i = 0; i < _seeds.size(); ++i)
    _weights[_seeds[i]] = i;
}

void Sector::prepare_targets(const std::vector<RawIntegral> &targets) {
  _generate_seeds();
}

unsigned Sector::run_reduce(const std::vector<IBPProtoFF> &ibps) {
  return sector_reduction(ibps);
}

unsigned Sector::run_reduce_sym(const std::vector<IBPProto> &ibps) {
  return sector_reduction_sym(ibps);
}

unsigned Sector::run_reduce_ff(const std::vector<IBPProto> &ibps) {
  _generate_seeds();
  return sector_reduction_ff(ibps);
}

unsigned Sector::sector_reduction(const std::vector<IBPProtoFF> &ibps) {
  // generate the system
  for (const auto &seed : _seeds) {
    if (seed.depth() < _depth && seed.rank() < _rank) {
      for (const auto &ibp : ibps) {
        EquationFF equation;
        // generate the ibp equation
        for (const auto &item : ibp) {
          RawIntegral integral = seed + item.first;
          if (!_weights.contains(integral))
            continue;
          // check if coefficient is zero
          umod64 coeff = item.second.back();
          for (unsigned k = 0; k < _nprops; ++k)
            if (item.second[k] != 0)
              coeff += item.second[k] * umod64::from(seed[k]);
          if (coeff == 0)
            continue;
          equation.insert(_weights[integral], coeff);
        }
        if (equation.empty())
          continue;
        else
          equation.sort();

        _systemFF.emplace_back(std::move(equation));
        _systemFF.back().eqnum = _systemFF.size();
      }
    }
  }
  std::sort(_systemFF.begin(), _systemFF.end());

  // gauss elimination
  for (auto &equation : _systemFF) {
    while (!equation.empty() && _lineNumber.contains(equation[0])) {
      equation.eliminate(_gaussFF[_lineNumber[equation[0]]], 0);
    }
    if (equation.empty())
      continue;
    equation.normalize();

    for (unsigned i = 1; i < equation.size();) {
      if (_lineNumber.contains(equation[i]))
        equation.eliminate(_gaussFF[_lineNumber[equation[i]]], i);
      else
        ++i;
    }

    _lineNumber[equation.first_integral()] = _gaussFF.size();
    _gaussFF.emplace_back(std::move(equation));
  }

  for (unsigned i = 0; i < _seeds.size(); ++i) {
    if (_seeds[i].depth() < _depth && _seeds[i].rank() < _rank) {
      if (!_lineNumber.contains(i))
        std::cout << "      " << _seeds[i] << "  # " << _id << std::endl;
    }
  }

  _systemFF.clear();
  _gaussFF.clear();
  _seeds.clear();
  _weights.clear();

  return 1;
}

unsigned Sector::sector_reduction_sym(const std::vector<IBPProto> &ibps) {
  // generate the system
  for (const auto &seed : _seeds) {
    if (seed.depth() < _depth && seed.rank() < _rank) {
      for (const auto &ibp : ibps) {
        EquationSym equation;
        // generate the ibp equation
        for (const auto &item : ibp) {
          RawIntegral integral = seed + item.first;
          if (!_weights.contains(integral))
            continue;
          // check if coefficient is zero
          GiNaC::lst values;
          for (unsigned p = 0; p < _nprops; ++p)
            values.append(_symIndices[p] == seed[p]);
          GiNaC::ex coeff = item.second.subs(values).expand();
          std::cout << coeff << std::endl;
          if (coeff == 0)
            continue;
          equation.insert(_weights[integral], coeff);
        }
        if (equation.empty())
          continue;
        else
          equation.sort();

        _systemS1.emplace_back(std::move(equation));
        _systemS1.back().eqnum = _systemS1.size();
      }
    }
  }
  std::sort(_systemS1.begin(), _systemS1.end());

  std::string recordPath = "record_" + std::to_string(_id);
  std::ofstream record(recordPath);

  // gauss elimination
  for (auto &equation : _systemS1) {
    std::cout << equation.eqnum << std::endl;
    while (!equation.empty() && _lineNumber.contains(equation[0])) {
      std::cout << _gaussS[_lineNumber[equation[0]]].eqnum << " ";
      equation.eliminate(_gaussS[_lineNumber[equation[0]]], 0);
    }
    if (equation.empty()) {
      std::cout << "\n" << std::endl;

      continue;
    }
    equation.normalize();

    for (unsigned i = 1; i < equation.size();) {
      if (_lineNumber.contains(equation[i])) {
        std::cout << _gaussS[_lineNumber[equation[i]]].eqnum << " ";
        equation.eliminate(_gaussS[_lineNumber[equation[i]]], i);
      } else
        ++i;
    }

    _lineNumber[equation.first_integral()] = _gaussS.size();
    _gaussS.emplace_back(std::move(equation));

    std::cout << "\n" << std::endl;
  }

  for (unsigned i = 0; i < _seeds.size(); ++i) {
    if (_seeds[i].depth() < _depth && _seeds[i].rank() < _rank) {
      if (!_lineNumber.contains(i))
        std::cout << "      " << _seeds[i] << "  # " << _id << std::endl;
    }
  }

  std::string path = "result_" + std::to_string(_id);
  std::ofstream file(path);

  for (const auto &item : _lineNumber) {
    file << _seeds[item.first] << std::endl;
    unsigned num = _gaussS[item.second].size();
    if (num > 1) {
      for (unsigned i = 1; i < num - 1; ++i)
        file << "(" << -_gaussS[item.second].coeff(i) << ")*"
             << _seeds[_gaussS[item.second][i]] << "+";
      file << "(" << -_gaussS[item.second].coeff(num - 1) << ")*"
           << _seeds[_gaussS[item.second][num - 1]];
      file << "\n" << std::endl;
    } else
      file << "0\n" << std::endl;
  }

  _systemS1.clear();
  _gaussS.clear();
  _seeds.clear();
  _weights.clear();

  return 1;
}

unsigned Sector::sector_reduction_ff(const std::vector<IBPProto> &ibps) {
  // clock_t startTime = clock();

  // std::cout << "  sector: " << _id << std::endl;
  // // create a new graph
  // Session session;
  // unsigned graphId = session.new_graph();
  // Graph *graph = session.graph(graphId);

  // // set input vars
  // graph->set_input_vars(0);
  // std::vector<unsigned> inputnodes;

  // // prepare sparse solver data
  // typedef NumericSparseSolverData Data;
  // std::unique_ptr<NumericSparseSolver> algptr(new NumericSparseSolver());
  // std::unique_ptr<Data> data(new Data());
  // auto &alg = *algptr;

  // std::ofstream file("result_" + std::to_string(_id));

  // unsigned nseeds = 0;
  // unsigned nouter = 0;
  // unsigned nsystem1 = 0;
  // unsigned nsystem2 = 0;
  // std::cout << "begin generating system" << std::endl;
  // std::set<RawIntegral> subSectors;
  // // generate the system
  // for (const auto &seed : _seeds) {
  //   if (seed.depth() < _depth && seed.rank() < _rank) {
  //     ++nseeds;
  //     for (const auto &ibp : ibps) {
  //       EquationSym equation;
  //       // generate the ibp equation
  //       for (const auto &item : ibp) {
  //         RawIntegral integral = seed + item.first;
  //         if (!_weights.contains(integral)) {
  //           subSectors.insert(integral);
  //           continue;
  //         }
  //         // check if coefficient is zero
  //         GiNaC::lst values;
  //         for (unsigned p = 0; p < _nprops; ++p)
  //           values.append(_symIndices[p] == seed[p]);
  //         GiNaC::ex coeff = item.second.subs(values).expand();
  //         if (coeff == 0)
  //           continue;
  //         equation.insert(_weights[integral], coeff);
  //       }
  //       if (equation.empty())
  //         continue;
  //       else
  //         equation.sort();

  //       RawIntegral firstInt = _seeds[equation.first_integral()];
  //       _systemS1.emplace_back(std::move(equation));
  //       _systemS1.back().eqnum = _systemS1.size();
  //       if (firstInt.depth() < _depth && firstInt.rank() < _rank) {
  //         ++nsystem1;
  //       } else
  //         ++nsystem2;
  //     }
  //   } else
  //     ++nouter;
  // }

  // // needed vars
  // unsigned nvars = _seeds.size();
  // unsigned neededvars = nseeds;
  // file << "integrals: " << nvars << std::endl;
  // file << "seeds: " << nseeds << "  outer: " << nouter << std::endl;
  // file << "system1: " << nsystem1 << std::endl;
  // file << "system2: " << nsystem2 << std::endl;
  // unsigned *neededVars = new unsigned[neededvars];
  // for (unsigned i = 0; i < neededvars; ++i)
  //   neededVars[i] = nouter + i;

  // // get system data
  // alg.rinfo.resize(_systemS1.size());
  // alg.c.resize(_systemS1.size());
  // for (unsigned i = 0; i < _systemS1.size(); ++i) {
  //   {
  //     unsigned csize = _systemS1[i].size();
  //     int *crinfo = new int[csize];
  //     for (unsigned it = 0; it < _systemS1[i].size(); ++it) {
  //       crinfo[it] = nvars - _systemS1[i][it] - 1;
  //     }

  //     NumericSparseSolver::RowInfo &rinf = alg.rinfo[i];
  //     rinf.size = csize;
  //     rinf.cols.reset(new unsigned[csize]);
  //     alg.c[i].reset(new MPRational[csize]);
  //     std::copy(crinfo, crinfo + csize, rinf.cols.get());
  //     delete[] crinfo;

  //     for (unsigned it = 0; it < _systemS1[i].size(); ++it) {
  //       std::stringstream ss;
  //       ss << _systemS1[i].coeff(it);
  //       alg.c[i][it] = MPRational(ss.str().c_str());
  //     }
  //   }
  // }

  // alg.init(alg.c.size(), nvars, neededVars, neededvars, *data);
  // delete[] neededVars;

  // // attach solver node
  // unsigned solverId =
  //     graph->new_node(std::move(algptr), std::move(data), inputnodes.data());

  // // set only homogeneous
  // Algorithm *presol = session.algorithm(graphId, solverId);
  // SparseLinearSolver &psls = *static_cast<SparseLinearSolver *>(presol);
  // psls.only_homogeneous();
  // session.invalidate_subctxt_alg_data(graphId, solverId);

  // // set output node
  // session.set_output_node(graphId, solverId);

  // clock_t genTime = clock();
  // file << "generate system in " << double(genTime - startTime) /
  // CLOCKS_PER_SEC
  //      << " s\n"
  //      << std::endl;

  // std::cout << "begin learn" << std::endl;
  // // solver learn
  // Ret lret = session.learn(graphId);
  // if (lret != SUCCESS) {
  //   std::cout << "learn failed!" << std::endl;
  //   exit(1);
  // }
  // // std::cout << "test" << std::endl;
  // clock_t learnTime = clock();
  // file << "find masters in " << double(learnTime - genTime) / CLOCKS_PER_SEC
  //      << " s" << std::endl;

  // // mark and sweep
  // Algorithm *solver = session.algorithm(graphId, solverId);
  // SparseLinearSolver &sls = *static_cast<SparseLinearSolver *>(solver);

  // file << "masters: ";
  // unsigned nmaster = 0;
  // for (unsigned i = 0; i < sls.n_needed_indepvars(); ++i) {
  //   RawIntegral integral = _seeds[nvars - 1 - sls.needed_indepvars()[i]];
  //   if (integral.depth() < _depth && integral.rank() < _rank) {
  //     nmaster++;
  //     file << _seeds[nvars - sls.needed_indepvars()[i] - 1] << " ";
  //   }
  // }
  // file << std::endl;
  // if (nmaster == 0) {
  //   std::cout << "subSector integrals: " << subSectors.size() << std::endl;
  //   _systemS1.clear();
  //   _systemS2.clear();
  //   _seeds.clear();
  //   _weights.clear();
  //   return 0;
  // }

  // sls.mark_and_sweep_eqs(session.alg_data(graphId, solverId));
  // session.invalidate_subctxt_alg_data(graphId, solverId);

  // unsigned nout = graph->nparsout;
  // file << "\nequations: " << sls.n_indep_eqs() << std::endl;
  // file << "reconstruct: " << nout << std::endl;

  // // delete unneeded equations
  // NumericSparseSolver &nss = *static_cast<NumericSparseSolver *>(solver);
  // nss.delete_unneeded_eqs(session.alg_data(graphId, solverId));
  // session.invalidate_subctxt_alg_data(graphId, solverId);

  // // reconstruction
  // ReconstructionOptions opt;
  // opt.max_primes = 100;

  // std::unique_ptr<MPRational[]> res(new MPRational[nout]);
  // Ret ret = session.reconstruct_numeric(graphId, res.get(), opt);
  // if (ret != SUCCESS) {
  //   std::cout << "reconstruction failed!" << std::endl;
  // }

  // clock_t reconsTime = clock();
  // file << "reconstruct in " << double(reconsTime - learnTime) /
  // CLOCKS_PER_SEC
  //      << " s" << std::endl;
  // file << "total time: " << double(reconsTime - startTime) / CLOCKS_PER_SEC
  //      << " s\n"
  //      << std::endl;

  // std::cout << "subSector integrals: " << subSectors.size() << "\n"
  //           << std::endl;

  // // result
  // for (unsigned i = 0; i < sls.n_needed_depvars(); ++i) {
  //   RawIntegral integral = _seeds[nvars - 1 - sls.needed_depvars()[i]];

  //   if (integral.depth() < _depth && integral.rank() < _rank) {
  //     file << "" << _seeds[nvars - 1 - sls.needed_depvars()[i]] << ":\n\n ";
  //     for (unsigned it = 0; it < sls.n_needed_indepvars(); ++it) {
  //       MPRational coeff = res[i * sls.n_needed_indepvars() + it];
  //       if (coeff != MPRational("0")) {
  //         file << res[i * sls.n_needed_indepvars() + it] << "*"
  //              << _seeds[nvars - 1 - sls.needed_indepvars()[it]] << "  ";
  //       }
  //     }
  //     file << "\n\n" << std::endl;
  //   }
  // }

  _systemS1.clear();
  _systemS2.clear();
  _seeds.clear();
  _weights.clear();
  _subSectors.clear();

  return 0;
}
