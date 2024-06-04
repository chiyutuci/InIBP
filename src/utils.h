#pragma once

#include <cstdlib>
#include <iostream>
#include <mutex>
#include <new>
#include <unordered_map>

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &p);

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  os << "[";
  for (auto it = v.begin(); it != v.end(); ++it) {
    if (it != v.begin())
      os << ", ";
    os << *it;
  }
  os << "]";
  return os;
}

template <typename T1, typename T2>
inline std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &p) {
  os << "<" << p.first << ", " << p.second << ">";
  return os;
}

inline void process_finish(const std::string &result) {
  std::cout << "\n\033[1m\033[36m Finish: \033[0m" << result << "\n"
            << std::endl;
  exit(EXIT_SUCCESS);
}
