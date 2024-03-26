#include "include/analysis.hpp"
#include "include/derivation.hpp"
#include "include/probability.hpp"
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <optional>
#include <set>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

const std::vector<Derivation> &Analysis::get_derivations() const {
  return derivations;
}

std::optional<Probability>
Analysis::get_relation_probability(const Relation &relation) const {
  auto it = relation_map.find(relation);
  if (it != relation_map.end()) {
    return it->second;
  }
  return std::nullopt;
}

std::optional<Probability>
Analysis::get_rule_probability(const Relation &head,
                               const std::vector<Relation> &body) const {
  std::vector<std::string> second;
  for (const auto &relation : body) {
    second.push_back(relation.name);
  }
  auto it = rule_map.find(Rule(head.name, second));
  if (it != rule_map.end()) {
    return it->second;
  }
  return std::nullopt;
}

void Analysis::calculate_probability() {
  bool fixed = false;
  std::map<std::pair<Relation, std::vector<Relation>>, std::vector<size_t>>
      rule_unknown_map{};
  std::map<std::pair<Relation, std::vector<Relation>>, Probability>
      rule_probability_map{};
  std::vector<Derivation> worklist;
  for (const auto &derivation : derivations) {
    if (relation_map.count(derivation.head) == 0) {
      worklist.push_back(derivation);
    }
  }
  while (!fixed) {
    fixed = true;
    for (const auto &derivation : worklist) {
      if (relation_map.count(derivation.head)) {
        continue;
      }
      std::optional<Probability> probability = Probability(0.0, 0.0);
      for (const auto &body : derivation.bodies) {
        auto rule_probability_it =
            rule_probability_map.find(std::make_pair(derivation.head, body));
        if (rule_probability_it != rule_probability_map.end()) {
          if (probability.has_value()) {
            probability = probability.value() | rule_probability_it->second;
          }
          continue;
        }

        Probability body_probability(1.0, 1.0);
        std::vector<size_t> unknown;

        for (size_t i = 0; i < body.size(); ++i) {
          auto relation_probability_it = relation_map.find(body[i]);
          if (relation_probability_it != relation_map.end()) {
            body_probability =
                body_probability & relation_probability_it->second;
          } else {
            unknown.push_back(i);
          }
        }

        if (unknown.empty()) {
          auto rule_probability = get_rule_probability(derivation.head, body);
          if (rule_probability.has_value()) {
            body_probability = body_probability * rule_probability.value();
          }
          if (probability.has_value()) {
            probability = probability.value() | body_probability;
          }
          rule_probability_map[std::make_pair(derivation.head, body)] =
              body_probability;
          fixed = false;
        } else {
          probability = std::nullopt;
          auto rule_unknown_it =
              rule_unknown_map.find(std::make_pair(derivation.head, body));
          if (rule_unknown_it != rule_unknown_map.end() &&
              rule_unknown_it->second == unknown) {
            continue;
          }
          fixed = false;
          rule_unknown_map[std::make_pair(derivation.head, body)] = unknown;
        }
      }
      if (probability.has_value()) {
        relation_map[derivation.head] = probability.value();
      }
    }
  }
}

void Analysis::dump(const std::optional<std::string> &path) const {
  if (path.has_value()) {
    std::ofstream ofs(path.value());
    if (!ofs.is_open()) {
      throw std::runtime_error("Failed to open file");
    } else {
      for (auto &[relation, probability] : relation_map) {
        ofs << relation.to_string() << ": " << probability.to_string()
            << std::endl;
      }
    }
    ofs.close();
  } else {
    for (auto &[relation, probability] : relation_map) {
      std::cout << relation.to_string() << ": " << probability.to_string()
                << std::endl;
    }
  }
}

Analysis::Analysis(const std::string &derivation_path,
                   const std::string &probability_path) {
  this->derivations = Derivation::load(derivation_path);
  std::tie(this->relation_map, this->rule_map) =
      Probability::load(probability_path);
}