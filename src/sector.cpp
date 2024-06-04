#include "sector.h"

#include <fstream>

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
          if (coeff == 0)
            continue;
          equation.insert(_weights[integral], coeff);
        }
        if (equation.empty())
          continue;
        else
          equation.sort();

        _systemS.emplace_back(std::move(equation));
        _systemS.back().eqnum = _systemS.size();
      }
    }
  }
  std::sort(_systemS.begin(), _systemS.end());

  std::string recordPath = "record_" + std::to_string(_id);
  std::ofstream record(recordPath);

  // gauss elimination
  for (auto &equation : _systemS) {
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

  _systemS.clear();
  _gaussS.clear();
  _seeds.clear();
  _weights.clear();

  return 1;
}