#include "include/probability.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

Probability Probability::operator&(const Probability &other) const {
  auto lower_bound = std::max(this->lower_bound + other.lower_bound - 1.0, 0.0);
  auto upper_bound = std::min(this->upper_bound, other.upper_bound);
  return Probability(lower_bound, upper_bound);
}

Probability Probability::operator|(const Probability &other) const {
  auto lower_bound = std::max(this->lower_bound, other.lower_bound);
  auto upper_bound = std::min(this->upper_bound + other.upper_bound, 1.0);
  return Probability(lower_bound, upper_bound);
}

Probability Probability::operator*(const Probability &other) const {
  auto lower_bound = this->lower_bound * other.lower_bound;
  auto upper_bound = this->upper_bound * other.upper_bound;
  return Probability(lower_bound, upper_bound);
}

std::tuple<std::map<Relation, Probability>, std::map<Rule, Probability>>
Probability::load(std::string path) {
  std::ifstream ifs(path);
  if (!ifs.is_open()) {
    throw std::runtime_error("Failed to open file: " + path);
  }

  std::map<Relation, Probability> relation_map;
  std::map<Rule, Probability> rule_map;

  std::string line;
  while (std::getline(ifs, line)) {
    auto terms = std::vector<std::string>();
    auto attrs = std::vector<std::string>();
    Probability prob;
    {
      std::istringstream iss(line);
      std::string term;
      while (iss >> term) {
        terms.push_back(term);
      }
      if (terms.size() != 4) {
        throw std::runtime_error("Invalid input file format");
      }
    }
    {
      std::istringstream iss(terms[2]);
      std::string attr;
      while (std::getline(iss, attr, ',')) {
        attrs.push_back(attr);
      }
    }
    if (terms[3].find(',') != std::string::npos) {
      std::istringstream iss(terms[3]);
      std::string lower_bound;
      std::string upper_bound;
      std::getline(iss, lower_bound, ',');
      std::getline(iss, upper_bound, ',');
      prob = Probability(std::stod(lower_bound), std::stod(upper_bound));
    } else {
      double p = std::stod(terms[3]);
      prob = Probability(p, p);
    }
    if (terms[0] == "relation") {
      std::vector<uint32_t> attrs_u32;
      for (const auto &attr : attrs) {
        attrs_u32.push_back(std::stoul(attr));
      }
      relation_map[Relation(terms[1], attrs_u32)] = prob;
    } else if (terms[0] == "rule") {
      rule_map[Rule(terms[1], attrs)] = prob;
    } else {
      throw std::runtime_error("Invalid input file format");
    }
  }
  ifs.close();
  return std::make_pair(relation_map, rule_map);
}

std::string Probability::to_string() const {
  std::stringstream ss;
  ss << "(" << std::fixed << std::setprecision(2) << lower_bound << ", " << upper_bound << ")";
  return ss.str();
}