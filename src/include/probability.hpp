#pragma once
#include "derivation.hpp"
#include <map>
#include <string>

struct Probability {
  double lower_bound = 0.0;
  double upper_bound = 1.0;

  Probability(double lower_bound = 0.0, double upper_bound = 1.0)
      : lower_bound(lower_bound), upper_bound(upper_bound) {}

  Probability operator&(const Probability &other) const;
  Probability operator|(const Probability &other) const;
  Probability operator*(const Probability &other) const;

  static std::tuple<std::map<Relation, Probability>, std::map<Rule, Probability>>
  load(std::string path);

  std::string to_string() const;
};
