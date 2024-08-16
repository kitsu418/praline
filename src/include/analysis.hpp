#pragma once
#include "derivation.hpp"
#include "probability.hpp"
#include <map>
#include <optional>
#include <vector>
#include <set>

class Analysis {
public:
  void calculate_probability(bool is_refined);
  void calculate_probability_legacy(bool is_refined);
  void dump(const std::optional<std::string> &path) const;
  void print_statistics() const;

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
  std::optional<Probability>
  get_relation_probability(const Relation &relation) const;
  std::optional<Probability>
  get_rule_probability(const Relation &head,
                       const std::vector<Relation> &body) const;
  void solve_unknowns();
  void refine();
  void strengthen_results();

private:
  std::vector<Derivation> derivations{};
  std::map<Relation, Probability> relation_map{};
  std::map<Rule, Probability> rule_map{};
  std::map<Relation, std::vector<Derivation>::iterator> derivations_index{};
  std::map<Relation, depclassid_t> dependent_class_map;
  std::map<std::pair<Relation, Relation>, Probability> dependent_map;
  std::set<std::string> queries;
  std::set<Relation> facts;
  bool is_legacy{false};

  std::set<Relation> dfs(Derivation &derivation, std::set<Derivation *> &visited,
           std::vector<Derivation *> &component);
};