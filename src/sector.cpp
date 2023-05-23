#include "sector.h"

std::map<std::pair<int, int>, std::vector<std::vector<int>>> Sector::combinations;

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
          for (const auto &comb: combinations[{num - 1, s - i}]) {
            std::vector<int> newComb = comb;
            newComb.emplace_back(i);
            combinations[{num, s}].emplace_back(std::move(newComb));
          }
        }
    }
}

void Sector::_generate_seeds() {
  int lines = std::popcount(_id);
  int zeros = (int) _nprops - lines;
  int number = std::max(lines, zeros);
  int sum = std::max((int) _depth - std::popcount(_id), (int) _rank);
  generate_combinations(number, sum);

  for (int depth = 0; depth <= _depth - std::popcount(_id); ++depth) {
    if (zeros != 0)
      for (int rank = 0; rank <= _rank; ++rank)
        for (const auto &depthComb: combinations[{lines, depth}])
          for (const auto &rankComb: combinations[{zeros, rank}]) {
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
      for (const auto &depthComb: combinations[{lines, depth}]) {
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

  for (const auto &target: targets)
    if (target.sector() == _id)
      _targets.insert(_weights[target]);

  for (auto target: _targets)
    _auxTargets.push(target);
  _usedTargets = _targets;
}

void Sector::run_reduce(const std::vector<IBPProto> &ibps,
                        const std::vector<GiNaC::symbol> &symbols) {
  // generate equations for each aux-target
  while (!_auxTargets.empty())
    _generate_equation(_auxTargets.front(), ibps, symbols);

  // find the master integrals
  std::cout << "\n number of master integrals: "
            << _usedTargets.size() - _system.size() << std::endl;
}

void Sector::_generate_equation(unsigned int target, const std::vector<IBPProto> &ibps,
                                const std::vector<GiNaC::symbol> &symbols) {
  RawIntegral integral(_seeds[target]);
  std::cout << "\n target: " << integral << std::endl;

  for (unsigned i = 0; i < ibps.size(); ++i) {
    for (const auto &item: ibps[i]) {
      // check if the seed is in this sector
      RawIntegral seed = integral - item.first;
      if (seed.sector() != _id)
        continue;
      // check if this ibp equation has been used
      // if not used, mark it as used
      if (_usedIBP.contains({_weights[seed], i + 1}))
        continue;
      else
        _usedIBP.insert({_weights[seed], i + 1});

      // get the value of indices
      GiNaC::lst values;
      for (unsigned p = 0; p < _nprops; ++p)
        values.append(symbols[p] == seed[p]);
      // check if this item is non-zero
      if (item.second.subs(values).expand() == 0)
        continue;
      // compute the equation
      std::map<unsigned, GiNaC::ex, std::greater<>> equation;
      bool allInSeeds = true;
      for (const auto &it: ibps[i]) {
        GiNaC::ex coeff = it.second.subs(values).expand();
        RawIntegral integ = seed + it.first;
        if (coeff != 0) {
          // check if the integral belongs to subsectors
          if (integ.sector() != _id)
            continue;
          // check if the integral belongs to _seeds
          if (_weights.contains(integ))
            equation.insert({_weights[integ], coeff});
          else {
            allInSeeds = false;
            break;
          }
        }
      }
      if (!allInSeeds || equation.empty())
        continue;

      // guass elimination
      while (!equation.empty() && _lineNumber.contains(equation.begin()->first)) {
        unsigned pivot = equation.begin()->first;
        GiNaC::ex scale = equation.begin()->second;
        for (const auto &term: _system[_lineNumber[pivot]])
          equation[term.first] -= scale * term.second;
        std::erase_if(equation, [](const auto &term) { return term.second.expand() == 0; });
      }
      if (equation.empty())
        continue;

      GiNaC::ex scale = equation.begin()->second;
      for (auto &term: equation)
        term.second /= scale;

      // record auxiliary targets
      for (const auto &term: equation)
        if (!_usedTargets.contains(term.first)) {
          _auxTargets.push(term.first);
          _usedTargets.insert(term.first);
        }

      std::cout << "  seed: " << seed << ", ibp: " << i + 1 << std::endl << "    ";
      for (const auto &term: equation)
        std::cout << term.second << " " << _seeds[term.first] << ", ";
      std::cout << std::endl;

      // add the equation to the system
      _lineNumber[equation.begin()->first] = _system.size();
      _system.emplace_back(std::move(equation));
      goto NextTarget;
    }
  }
  NextTarget:
  _auxTargets.pop();
}
