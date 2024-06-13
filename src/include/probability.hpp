#pragma once
#include "derivation.hpp"
#include <map>
#include <set>
#include <string>

struct Probability {
  double lower_bound = 0.0;
  double upper_bound = 1.0;

  Probability(double lower_bound = 0.0, double upper_bound = 1.0)
      : lower_bound(lower_bound), upper_bound(upper_bound) {}

  Probability dep_conj(const Probability &other) const;
  Probability dep_disj(const Probability &other) const;
  Probability pos_conj(const Probability &other) const;
  Probability pos_disj(const Probability &other) const;
  Probability neg_conj(const Probability &other) const;
  Probability neg_disj(const Probability &other) const;
  Probability ind_conj(const Probability &other) const;
  Probability ind_disj(const Probability &other) const;
  
  static std::tuple<std::map<Relation, Probability>,
                    std::map<Rule, Probability>,
                    std::map<Relation, depclassid_t>,
                    std::map<std::pair<Relation, Relation>, Probability>,
                    std::set<std::string>, std::set<Relation>>
  load(const std::string &path, const bool &is_legacy);

  std::string to_string() const;
};
