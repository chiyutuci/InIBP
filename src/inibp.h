#pragma once

#include <iostream>

#include <yaml-cpp/yaml.h>

#include "family.h"

class InIBP {
public:
  explicit InIBP(const YAML::Node &);

public:
  // run the reduction
  void run();

private:
  Family _family;
  Reduce _reduce;
};