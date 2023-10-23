#include "sector.h"
#include <ctime>

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

unsigned Sector::block_reduce(unsigned int depth, unsigned int rank,
                          const std::vector<IBPProtoFF> &ibps) {
  std::vector<RawIntegral> seeds;
  std::vector<EquationFF> system;
  std::unordered_map<unsigned,unsigned> lineNumber;

  // generate all seeds within this block
  int lines = std::popcount(_id);
  int zeros = (int) _nprops - lines;
  if (zeros!=0){
    for (const auto &depthComb: combinations[{lines, depth-lines}])
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
        seeds.emplace_back(std::move(integral));
      }
  }
  else {
    for (const auto &depthComb: combinations[{lines, depth-lines}]) {
      RawIntegral integral(_lines);
      unsigned posLines = 0;
      for (unsigned i = 0; i < _nprops; i++)
        if (_lines[i])
          integral[i] += depthComb[posLines++];
      seeds.emplace_back(std::move(integral));
    }
  }

  for (const auto& seed: seeds) {
    for (const auto & ibp : ibps) {
      EquationFF equation;
      //generate the ibp equation
      for (const auto& item: ibp) {
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

      // Gauss elimination
      while (!equation.empty() && lineNumber.contains(equation.first_integral())) {
        equation.eliminate(system[lineNumber[equation.first_integral()]], 0);
      }
      if (equation.empty())
        continue;
      equation.normalize();
      // add this equation to the system
      lineNumber[equation.first_integral()] = system.size();
      system.emplace_back(std::move(equation));
    }
  }

  return system.size();
}

void Sector::prepare_targets(const std::vector<RawIntegral> &targets) {
  _generate_seeds();
}

unsigned Sector::run_reduce(const std::vector<IBPProtoFF> &ibps) {
//  std::cout << "\n  sector: " << _id << "\n" << std::endl;
//  unsigned equations = 0;
//    for (unsigned depth = std::popcount(_id); depth < _depth; ++depth)
//      for (unsigned rank = 0; rank < _rank; ++rank)
//        equations += block_reduce(depth, rank, ibps);
//  return equations;
  return sector_reduction(ibps);
}

unsigned Sector::sector_reduction(const std::vector<IBPProtoFF> & ibps) {
  // generate the system
  for (const auto& seed: _seeds) {
    if (seed.depth() < _depth && seed.rank() < _rank) {
      for (const auto &ibp: ibps) {
        EquationFF equation;
        //generate the ibp equation
        for (const auto &item: ibp) {
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
  for (auto& equation : _systemFF) {
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
        std::cout << "      " << _seeds[i] << std::endl;
    }
  }

  return _gaussFF.size();
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
        equation.eliminate(_systemFF[_lineNumber[equation.first_integral()]],0);
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
