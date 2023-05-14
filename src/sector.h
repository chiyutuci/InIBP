#pragma once


class Sector {
public:
  [[nodiscard]] unsigned id() const { return _id; }

  unsigned &id() { return _id; }

private:
  unsigned _id;
};