#pragma once

#include <vector>

class Reduce;

class Sector {
public:
  friend class Reduce;

private:
  // sector number
  unsigned _id;
  // super sectors
  std::vector<unsigned> _superSectors;
  // sub sectors
  std::vector<unsigned> _subSectors;
};