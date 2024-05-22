#pragma once
#include <string>
#include <vector>

using depclassid_t = uint32_t;

struct Relation {
  std::string name;
  std::vector<uint32_t> attributes;

  Relation(std::string name, std::vector<uint32_t> attributes)
      : name(name), attributes(attributes) {}
  Relation() = default;

  bool operator<(const Relation &other) const;
  bool operator==(const Relation &other) const;
  bool operator!=(const Relation &other) const;

  std::string to_string() const;
};

struct Rule {
  std::string head_name;
  std::vector<std::string> body_name;

  Rule(std::string head_name, std::vector<std::string> body_name)
      : head_name(head_name), body_name(body_name) {}
  Rule() = default;

  bool operator<(const Rule &other) const;
  bool operator==(const Rule &other) const;
};

struct Derivation {
  Relation head;
  std::vector<std::vector<Relation>> bodies;

  Derivation(Relation head, std::vector<std::vector<Relation>> bodies)
      : head(head), bodies(bodies) {}
  Derivation() = default;

  static std::vector<Derivation> load(const std::string &path);
  std::string to_string() const;

  bool operator<(const Derivation &other) const;
};
