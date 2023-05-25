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

void Sector::run_reduce(const std::vector<IBPProtoFF> &ibps) {
  // generate equations for each aux-target
  while (!_auxTargets.empty())
    _generate_equation(_auxTargets.front(), ibps);

  // find the master integrals
  unsigned nmis = 0;
  std::cout << "\n\n master integrals:\n" << std::endl;
  for (unsigned integral = 0; integral <= *_targets.begin(); ++integral)
    if (!_lineNumber.contains(integral)) {
      std::cout << "   " << _seeds[integral] << std::endl;
      ++nmis;
    }
  std::cout << "number of master integrals: " << nmis << std::endl;

  unsigned neqs = 0;
  for (const auto &item: _lineNumber) {
    if (item.first > *_targets.begin())
      ++neqs;
  }
  std::cout << "\n number of equations: " << _systemFF.size() << std::endl;
  std::cout << " out bound: " << neqs << std::endl;

}

void Sector::_generate_equation(unsigned int target, const std::vector<IBPProtoFF> &ibps) {
  RawIntegral integral(_seeds[target]);

  for (unsigned i = 0; i < ibps.size(); ++i) {
    for (const auto &item: ibps[i]) {
      // check if the seed is in this sector
      RawIntegral seed = integral - item.first;
      if (seed.sector() != _id || !_weights.contains(seed))
        continue;
      // check if this ibp equation has been used
      // if not used, mark it as used
      if (_usedIBP.contains({_weights[seed], i + 1}))
        continue;

      // check if this item is non-zero
      umod64 coe = item.second.back();
      for (unsigned k = 0; k < _nprops; ++k)
        if (item.second[k] != 0)
          coe += item.second[k] * umod64::from(seed[k]);
      if (coe == 0)
        continue;

      // compute the equation
      EquationFF equation;
      bool allInSeeds = true;
      for (const auto &it: ibps[i]) {
        RawIntegral integ = seed + it.first;
        // check if the integral belongs to subsectors
        if (integ.sector() != _id)
          continue;
        umod64 coeff = it.second.back();
        for (unsigned k = 0; k < _nprops; ++k)
          if (it.second[k] != 0)
            coeff += it.second[k] * umod64::from(seed[k]);
        if (coeff != 0) {
          // check if the integral belongs to _seeds
          if (_weights.contains(integ)) {
            equation.insert(_weights[integ], coeff);
          }
          else {
            allInSeeds = false;
            break;
          }
        }
      }
      if (!allInSeeds || equation.empty()) {
        _usedIBP.insert({_weights[seed], i + 1});
        continue;
      }

      // guass elimination
      while (!equation.empty() && _lineNumber.contains(equation.first_integral())) {
        equation.eliminate(_systemFF[_lineNumber[equation.first_integral()]]);
      }
      if (equation.empty()) {
        _usedIBP.insert({_weights[seed], i + 1});
        continue;
      }

      equation.normalize();

      // record auxiliary targets
      for (const auto &term: equation.eq())
        if (!_usedTargets.contains(term.first)) {
          _auxTargets.push(term.first);
          _usedTargets.insert(term.first);
        }

      // add the equation to the system
      _lineNumber[equation.first_integral()] = _systemFF.size();
      _systemFF.emplace_back(std::move(equation));
      _usedIBP.insert({_weights[seed], i + 1});
      // goto NextTarget;
    }
  }
  NextTarget:
  _auxTargets.pop();
}
