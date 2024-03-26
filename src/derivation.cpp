#include "include/derivation.hpp"
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

bool Relation::operator<(const Relation &other) const {
  if (name != other.name) {
    return name < other.name;
  }
  return attributes < other.attributes;
}

bool Relation::operator==(const Relation &other) const {
  return name == other.name && attributes == other.attributes;
}

bool Relation::operator!=(const Relation &other) const {
  return !(*this == other);
}

std::string Relation::to_string() const {
  std::string str = name + "(";
  for (size_t i = 0; i < attributes.size(); i++) {
    str += std::to_string(attributes[i]);
    if (i != attributes.size() - 1) {
      str += ", ";
    }
  }
  str += ")";
  return str;
}

bool Rule::operator<(const Rule &other) const {
  if (head_name != other.head_name) {
    return head_name < other.head_name;
  }
  return body_name < other.body_name;
}

bool Rule::operator==(const Rule &other) const {
  return head_name == other.head_name && body_name == other.body_name;
}

using json = nlohmann::json;

static void from_json(const json &j, Relation &r) {
  j.at("name").get_to(r.name);
  j.at("attributes").get_to(r.attributes);
}

static void from_json(const json &j, std::vector<Relation> &v) {
  v.clear();
  for (const auto &relation : j) {
    Relation r;
    relation.get_to(r);
    v.push_back(r);
  }
}

static void from_json(const json &j,
                      std::vector<std::vector<Relation>> &nested_vec) {
  nested_vec.clear();
  for (const auto &inner : j) {
    std::vector<Relation> inner_vec;
    inner.get_to(inner_vec);
    nested_vec.push_back(inner_vec);
  }
}

static void from_json(const json &j, Derivation &d) {
  j.at("head").get_to(d.head);
  j.at("bodies").get_to(d.bodies);
}

static void from_json(const json &j, std::vector<Derivation> &v) {
  v.clear();
  for (const auto &derivation : j) {
    Derivation d;
    derivation.get_to(d);
    v.push_back(d);
  }
}

std::vector<Derivation> Derivation::load(const std::string &path) {
  std::ifstream ifs(path);
  if (!ifs.is_open()) {
    throw std::runtime_error("Failed to open file: " + path);
  }

  json j;
  ifs >> j;
  ifs.close();

  std::vector<Derivation> derivations;
  from_json(j, derivations);

  return derivations;
}