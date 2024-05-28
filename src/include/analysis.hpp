#pragma once
#include "derivation.hpp"
#include "probability.hpp"
#include <map>
#include <optional>
#include <vector>

class Analysis {
public:
  void calculate_probability();
  void calculate_probability_legacy();
  void dump(const std::optional<std::string> &path) const;

public:
  Analysis(const std::string &derivation_path,
           const std::string &probability_path, const bool &is_legcy);

private:
  Analysis(const std::vector<Derivation> &derivations,
           const std::map<Relation, Probability> &relation_map,
           const std::map<Rule, Probability> &rule_map,
           const std::map<Relation, std::vector<Derivation>::iterator>
               &derivations_index)
      : derivations(derivations), relation_map(relation_map),
        rule_map(rule_map), derivations_index(derivations_index){};
  const std::vector<Derivation> &get_derivations() const;
  std::optional<Probability>
  get_relation_probability(const Relation &relation) const;
  std::optional<Probability>
  get_rule_probability(const Relation &head,
                       const std::vector<Relation> &body) const;
  void solve_unknowns();
  void compute(const Relation &head);

private:
  std::vector<Derivation> derivations{};
  std::map<Relation, Probability> relation_map{};
  std::map<Rule, Probability> rule_map{};
  std::map<Relation, std::vector<Derivation>::iterator> derivations_index{};
  std::map<Relation, depclassid_t> dependent_class_map;
  std::map<std::pair<Relation, Relation>, Probability> dependent_map;
  bool is_legacy{false};
};