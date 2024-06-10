#include "include/probability.hpp"
#include "include/derivation.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

Probability Probability::dep_conj(const Probability &other) const {
  auto lower_bound = std::max(this->lower_bound + other.lower_bound - 1.0, 0.0);
  auto upper_bound = std::min(this->upper_bound, other.upper_bound);
  return Probability(std::move(lower_bound), std::move(upper_bound));
}

Probability Probability::dep_disj(const Probability &other) const {
  auto lower_bound = std::max(this->lower_bound, other.lower_bound);
  auto upper_bound = std::min(this->upper_bound + other.upper_bound, 1.0);
  return Probability(std::move(lower_bound), std::move(upper_bound));
}

Probability Probability::pos_conj(const Probability &other) const {
  auto lower_bound = this->lower_bound * other.lower_bound;
  auto upper_bound = std::min(this->upper_bound, other.upper_bound);
  return Probability(std::move(lower_bound), std::move(upper_bound));
}

Probability Probability::pos_disj(const Probability &other) const {
  auto lower_bound = std::max(this->lower_bound, other.lower_bound);
  auto upper_bound = 1.0 - (1 - this->upper_bound) * (1 - other.upper_bound);
  return Probability(std::move(lower_bound), std::move(upper_bound));
}

Probability Probability::neg_conj(const Probability &other) const {
  auto lower_bound = std::max(this->lower_bound + other.lower_bound - 1.0, 0.0);
  auto upper_bound = this->upper_bound * other.upper_bound;
  return Probability(std::move(lower_bound), std::move(upper_bound));
}

Probability Probability::neg_disj(const Probability &other) const {
  auto lower_bound =
      1.0 - (1.0 - this->lower_bound) * (1.0 - other.lower_bound);
  auto upper_bound = std::min(1.0, this->upper_bound + other.upper_bound);
  return Probability(std::move(lower_bound), std::move(upper_bound));
}

Probability Probability::ind_conj(const Probability &other) const {
  auto lower_bound = this->lower_bound * other.lower_bound;
  auto upper_bound = this->upper_bound * other.upper_bound;
  return Probability(std::move(lower_bound), std::move(upper_bound));
}

Probability Probability::ind_disj(const Probability &other) const {
  auto lower_bound =
      1.0 - (1.0 - this->lower_bound) * (1.0 - other.lower_bound);
  auto upper_bound =
      1.0 - (1.0 - this->upper_bound) * (1.0 - other.upper_bound);
  return Probability(std::move(lower_bound), std::move(upper_bound));
}

std::tuple<std::map<Relation, Probability>, std::map<Rule, Probability>,
           std::map<Relation, depclassid_t>,
           std::map<std::pair<Relation, Relation>, Probability>,
           std::set<std::string>>
Probability::load(const std::string &path, const bool &is_legacy) {
  std::ifstream ifs(path);
  if (!ifs.is_open()) {
    throw std::runtime_error("Failed to open file: " + path);
  }

  std::map<Relation, Probability> relation_map;
  std::map<Rule, Probability> rule_map;
  std::map<Relation, depclassid_t> dependent_class_map;
  std::map<std::pair<Relation, Relation>, Probability> dependent_map;
  std::set<std::string> queries;

  auto parse_attributes = []<typename T>(const std::string &line,
                                         std::vector<T> &attrs) {
    std::istringstream iss(line);
    std::string attr;
    while (std::getline(iss, attr, ',')) {
      if constexpr (std::is_same<T, uint32_t>::value) {
        attrs.push_back(std::stoul(attr));
      } else if (std::is_same<T, std::string>::value) {
        attrs.push_back(attr);
      }
    }
    return attrs;
  };

  auto parse_probability = [](const std::string &line) {
    if (line.find(',') != std::string::npos) {
      std::istringstream iss(line);
      std::string lower_bound;
      std::string upper_bound;
      std::getline(iss, lower_bound, ',');
      std::getline(iss, upper_bound, ',');
      return Probability(std::stod(lower_bound), std::stod(upper_bound));
    } else {
      double p = std::stod(line);
      return std::move(Probability(p, p));
    }
  };

  std::string line;
  while (std::getline(ifs, line)) {
    // relation: 0: type, 1: name, 2: attributes, 3: probability, 4: class
    // rule: 0: type, 1: name, 2: attributes, 3: probability
    // dependent: 0: type, 1: relation_name1, 2: attributes_1, 3:
    // relation_name2, 4: attributes_2, 5: probability
    auto terms = std::vector<std::string>();
    Probability prob;

    {
      std::stringstream iss(line);
      std::string term;
      while (std::getline(iss, term, ';')) {
        terms.push_back(term);
      }
    }

    if (terms[0] == "rule" || terms[0] == "relation") {
      if (terms[0] == "relation") {
        if ((terms.size() != 5 && !is_legacy) ||
            (is_legacy && terms.size() != 4 && terms.size() != 5)) {
          throw std::runtime_error("Invalid input file format (relation): " +
                                   line);
        }
        std::vector<uint32_t> attrs;
        parse_attributes(terms[2], attrs);
        auto relation = Relation(terms[1], std::move(attrs));
        relation_map[relation] = std::move(parse_probability(terms[3]));
        if (!is_legacy) {
          auto depclass_id = std::stoul(terms[4]);
          dependent_class_map[relation] = depclass_id;
        }
      } else {
        if (terms.size() != 4) {
          throw std::runtime_error("Invalid input file format (rule): " + line);
        }
        std::vector<std::string> attrs;
        parse_attributes(terms[2], attrs);
        rule_map[Rule(terms[1], std::move(attrs))] =
            std::move(parse_probability(terms[3]));
      }
    } else if (terms[0] == "query") {
      if (terms.size() != 2) {
        throw std::runtime_error("Invalid input file format (query): " + line);
      }
      queries.insert(terms[1]);
    } else if (!is_legacy && terms[0] == "dependent") {
      if (is_legacy || terms.size() != 6) {
        throw std::runtime_error("Invalid input file format (dependent): " +
                                 line);
      }
      std::vector<uint32_t> attrs1;
      parse_attributes(terms[2], attrs1);
      std::vector<uint32_t> attrs2;
      parse_attributes(terms[4], attrs2);
      dependent_map[std::make_pair(Relation(terms[1], std::move(attrs1)),
                                   Relation(terms[3], std::move(attrs2)))] =
          std::stod(terms[5]);
    }
  }
  ifs.close();
  return std::move(std::make_tuple(std::move(relation_map), std::move(rule_map),
                                   std::move(dependent_class_map),
                                   std::move(dependent_map),
                                   std::move(queries)));
}

std::string Probability::to_string() const {
  std::stringstream ss;
  ss << "[" << std::fixed << std::setprecision(6) << lower_bound << ","
     << upper_bound << "]";
  return ss.str();
}