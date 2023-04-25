#pragma once

#include <iostream>
#include <vector>


template<typename T1, typename T2>
inline std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &p);

template<typename T>
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

template<typename T1, typename T2>
inline std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &p) {
	os << "<" << p.first << ", " << p.second << ">";
	return os;
}