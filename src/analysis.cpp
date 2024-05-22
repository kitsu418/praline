#include "include/analysis.hpp"
#include "include/derivation.hpp"
#include "include/probability.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cstddef>
#include <cstdio>
#include <cvc5/cvc5.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <tuple>
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

enum Operation {
  kConjunction,
  kDisjunction,
};

// void Analysis::compute(const Relation &head) {
//   auto head_it = derivations_index[head];
//   Probability probability(0.0, 0.0);
//   for (const auto &body : head_it->bodies) {
//     Probability body_probability(1.0, 1.0);
//     for (const auto &relation : body) {
//       auto relation_probability_it = relation_map.find(relation);
//       if (relation_probability_it != relation_map.end()) {
//         body_probability =
//             body_probability.dep_conj(relation_probability_it->second);
//       } else {
//         compute(relation);
//       }
//     }
//     auto rule_probability = get_rule_probability(head, body);
//     if (rule_probability.has_value()) {
//       body_probability = body_probability.ind_conj(rule_probability.value());
//     }
//     probability = probability.dep_disj(body_probability);
//   }
//   relation_map[head] = std::move(probability);
// }

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

  std::map<Relation, std::set<Relation>> relation_dep_leaves;
  std::map<std::pair<Relation, Relation>, bool> subset_is_ind;

  std::function<std::set<Relation>(const Relation &,
                                   std::map<Relation, bool> &)>
      collect = [&](const Relation &relation,
                    std::map<Relation, bool> &visited) {
        std::set<Relation> leaves;
        if (dependent_class_map.find(relation) != dependent_class_map.end()) {
          leaves.insert(relation);
        } else {
          if (!visited[relation]) {
            visited[relation] = true;
            for (const auto &body : derivations_index[relation]->bodies) {
              for (const auto &r : body) {
                std::set<Relation> tmp = collect(r, visited);
                leaves.merge(tmp);
              }
            }
          }
        }
        return std::move(leaves);
      };

  auto operate = [&](const std::set<Relation> &lhs, const Probability &lhsp,
                     const std::set<Relation> &rhs, const Probability &rhsp,
                     const Operation &op) {
    std::set<Relation> lset;
    std::set<Relation> rset;

    for (auto &relation : lhs) {
      auto relation_dep_leaves_it = relation_dep_leaves.find(relation);
      if (relation_dep_leaves_it == relation_dep_leaves.end()) {
        std::map<Relation, bool> visited;
        auto collected = collect(relation, visited);
        relation_dep_leaves[relation] = collected;
        lset.merge(collected);
      } else {
        for (auto &r : relation_dep_leaves[relation]) {
          lset.insert(r);
        }
      }
    }
    for (auto &relation : rhs) {
      auto relation_dep_leaves_it = relation_dep_leaves.find(relation);
      if (relation_dep_leaves_it == relation_dep_leaves.end()) {
        std::map<Relation, bool> visited;
        auto collected = collect(relation, visited);
        relation_dep_leaves[relation] = collected;
        rset.merge(collected);
      } else {
        for (auto &r : relation_dep_leaves[relation]) {
          rset.insert(r);
        }
      }
    }

    std::set<Relation> shared;
    for (auto &r : lset) {
      if (rset.contains(r)) {
        shared.insert(r);
      }
    }

    std::set<depclassid_t> lids;
    std::set<depclassid_t> rids;

    for (auto &relation : lset) {
      if (!shared.contains(relation)) {
        lids.insert(dependent_class_map[relation]);
      }
    }
    for (auto &relation : rset) {
      if (!shared.contains(relation)) {
        rids.insert(dependent_class_map[relation]);
      }
    }

    for (auto &l1 : lids) {
      for (auto &l2 : rids) {
        if (l1 == l2) {
          if (op == kConjunction) {
            return lhsp.dep_conj(rhsp);
          } else {
            return lhsp.dep_disj(rhsp);
          }
        }
      }
    }

    if (shared.empty()) {
      if (op == kConjunction) {
        return lhsp.ind_conj(rhsp);
      } else {
        return lhsp.ind_disj(rhsp);
      }
    } else {
      if (op == kConjunction) {
        return lhsp.pos_conj(rhsp);
      } else {
        return lhsp.pos_disj(rhsp);
      }
    }
  };

  while (!fixed) {
    fixed = true;
    for (const auto &derivation : worklist) {
      if (relation_map.count(derivation.head)) {
        continue;
      }

      std::optional<Probability> probability = Probability(0.0, 0.0);
      std::set<Relation> lhsd;

      for (const auto &body : derivation.bodies) {
        auto rule_probability_it =
            rule_probability_map.find(std::make_pair(derivation.head, body));
        if (rule_probability_it != rule_probability_map.end()) {
          if (probability.has_value()) {
            probability =
                operate(lhsd, probability.value(), {body.begin(), body.end()},
                        rule_probability_it->second, kDisjunction);
            for (auto &r : body) {
              lhsd.insert(r);
            }
          }
          continue;
        }

        Probability body_probability(1.0, 1.0);
        std::vector<size_t> unknown;
        std::set<Relation> lhsc;

        for (size_t i = 0; i < body.size(); ++i) {
          auto relation_probability_it = relation_map.find(body[i]);
          if (relation_probability_it != relation_map.end() &&
              find(lhsc.begin(), lhsc.end(), body[i]) == lhsc.end()) {
            for (int j = i + 1; j < body.size(); ++j) {
              auto condp_it =
                  dependent_map.find(std::make_pair(body[i], body[j]));
              if (condp_it != dependent_map.end()) {
                body_probability = operate(
                    lhsc, body_probability, {body[i], body[j]},
                    condp_it->second.ind_conj(relation_probability_it->second),
                    kConjunction);
                lhsc.insert(body[j]);
                goto next_rel;
              }
            }
            body_probability =
                operate(lhsc, body_probability, {body[i]},
                        relation_probability_it->second, kConjunction);
          next_rel:
            lhsc.insert(body[i]);
          } else {
            unknown.push_back(i);
            lhsc.insert(body[i]);
          }
        }

        if (unknown.empty()) {
          auto rule_probability = get_rule_probability(derivation.head, body);
          if (rule_probability.has_value()) {
            body_probability =
                body_probability.ind_conj(rule_probability.value());
          }
          if (probability.has_value()) {
            probability = probability.value().dep_disj(body_probability);
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

        lhsd.merge(lhsc);
      }
      if (probability.has_value()) {
        relation_map[derivation.head] = probability.value();
      }
    }
  }
  solve_unknowns();
}

void Analysis::dump(const std::optional<std::string> &path) const {
  if (path.has_value()) {
    std::ofstream ofs(path.value());
    if (!ofs.is_open()) {
      throw std::runtime_error("Failed to open file");
    } else {
      for (auto &[relation, probability] : relation_map) {
        ofs << relation.to_string() << " " << probability.to_string()
            << std::endl;
      }
    }
    ofs.close();
  } else {
    for (auto &[relation, probability] : relation_map) {
      std::cout << relation.to_string() << " " << probability.to_string()
                << std::endl;
    }
  }
}

Analysis::Analysis(const std::string &derivation_path,
                   const std::string &probability_path) {
  this->derivations = Derivation::load(derivation_path);
  std::tie(this->relation_map, this->rule_map, this->dependent_class_map,
           this->dependent_map) = Probability::load(probability_path);
  for (auto it = this->derivations.begin(); it != this->derivations.end();
       ++it) {
    this->derivations_index[it->head] = it;
  }
}

using namespace cvc5;

void Analysis::solve_unknowns() {
  TermManager tm;
  Solver solver(tm);
  solver.setOption("produce-models", "true");
  solver.setOption("produce-unsat-cores", "true");
  solver.setLogic("ALL");

  Term rm = tm.mkRoundingMode(RoundingMode::ROUND_NEAREST_TIES_TO_EVEN);
  Sort fpt64 = tm.mkFloatingPointSort(11, 53);

  auto term_name_generator = [](const Relation &relation) {
    std::stringstream ss;
    ss << relation.name;
    for (const auto &attribute : relation.attributes) {
      ss << "_" << attribute;
    }
    return ss.str();
  };

  auto get_fpt64_value = [](const Term &term) -> double {
    uint64_t bits = std::stoull(
        std::get<2>(term.getFloatingPointValue()).getBitVectorValue(10));
    return std::bit_cast<double>(bits);
  };

  auto mk_fpt64_constant = [&](const double &value) {
    return tm.mkFloatingPoint(
        11, 53, tm.mkBitVector(64, std::bit_cast<uint64_t>(value)));
  };

  auto zero = mk_fpt64_constant(0.0);
  auto one = mk_fpt64_constant(1.0);
  std::map<Relation, std::array<Term, 2>> relation_terms;

  for (const auto &derivation : derivations) {
    if (relation_map.find(derivation.head) != relation_map.end()) {
      continue;
    };

    std::queue<Relation> worklist;
    std::set<Relation> asserted;
    worklist.push(derivation.head);
    solver.resetAssertions();
    asserted.insert(derivation.head);

    while (!worklist.empty()) {
      auto head = worklist.front();
      worklist.pop();

      if (relation_terms.find(head) == relation_terms.end()) {
        relation_terms[head] = std::to_array(
            {tm.mkConst(fpt64, term_name_generator(head) + "_l"),
             tm.mkConst(fpt64, term_name_generator(head) + "_u")});
      }

      Term head_term[2] = {relation_terms[head][0], relation_terms[head][1]};
      Term result_term[2];
      std::vector<Term> body_terms[2];
      for (const auto &body : derivations_index[head]->bodies) {
        Term body_term[2];
        std::vector<Term> terms[2];
        for (const auto &relation : body) {
          Term term[2];
          if (relation_map.find(relation) != relation_map.end()) {
            term[0] = mk_fpt64_constant(relation_map[relation].lower_bound);
            term[1] = mk_fpt64_constant(relation_map[relation].upper_bound);
          } else {
            if (!asserted.contains(relation)) {
              asserted.insert(relation);
              worklist.push(relation);
            }
            if (relation_terms.find(relation) == relation_terms.end()) {
              relation_terms[relation] = std::to_array(
                  {tm.mkConst(fpt64, term_name_generator(relation) + "_l"),
                   tm.mkConst(fpt64, term_name_generator(relation) + "_u")});
            }
            for (size_t i = 0; i < 2; ++i) {
              term[i] = relation_terms[relation][i];
            }
          }
          for (size_t i = 0; i < 2; ++i) {
            terms[i].push_back(term[i]);
          }
        }
        if (terms[0].size() < 1) {
          continue;
        }
        body_term[0] = terms[0][0];
        body_term[1] = terms[1][0];
        for (size_t i = 1; i < terms[0].size(); ++i) {
          auto tmp = tm.mkTerm(Kind::FLOATINGPOINT_SUB,
                               {rm,
                                tm.mkTerm(Kind::FLOATINGPOINT_ADD,
                                          {rm, body_term[0], terms[0][i]}),
                                one});
          body_term[0] = tm.mkTerm(
              Kind::ITE,
              {tm.mkTerm(Kind::FLOATINGPOINT_GT, {tmp, zero}), tmp, zero});
          body_term[1] = tm.mkTerm(
              Kind::ITE,
              {tm.mkTerm(Kind::FLOATINGPOINT_LT, {body_term[1], terms[1][i]}),
               body_term[1], terms[1][i]});
        }
        auto rule_p = get_rule_probability(head, body);
        if (rule_p.has_value()) {
          for (size_t i = 0; i < 2; ++i) {
            body_term[i] = tm.mkTerm(
                Kind::FLOATINGPOINT_MULT,
                {rm, body_term[i],
                 mk_fpt64_constant(i == 0 ? rule_p.value().lower_bound
                                          : rule_p.value().upper_bound)});
          }
        }
        for (size_t i = 0; i < 2; ++i) {
          body_terms[i].push_back(body_term[i]);
        }
      }
      if (body_terms[0].size() < 1) {
        continue;
      }
      result_term[0] = body_terms[0][0];
      result_term[1] = body_terms[1][0];
      for (size_t i = 1; i < body_terms[0].size(); ++i) {
        result_term[0] =
            tm.mkTerm(Kind::ITE, {tm.mkTerm(Kind::FLOATINGPOINT_GT,
                                            {result_term[0], body_terms[0][i]}),
                                  result_term[0], body_terms[0][i]});
        auto tmp = tm.mkTerm(Kind::FLOATINGPOINT_ADD,
                             {rm, result_term[1], body_terms[1][i]});
        result_term[1] =
            tm.mkTerm(Kind::ITE, {tm.mkTerm(Kind::FLOATINGPOINT_LT, {tmp, one}),
                                  tmp, one});
      }
      for (int i = 0; i < 2; ++i) {
        solver.assertFormula(
            tm.mkTerm(Kind::FLOATINGPOINT_EQ, {head_term[i], result_term[i]}));
        solver.assertFormula(
            tm.mkTerm(Kind::FLOATINGPOINT_LEQ, {head_term[i], one}));
        solver.assertFormula(
            tm.mkTerm(Kind::FLOATINGPOINT_GEQ, {head_term[i], zero}));
      }
    }
    Result r = solver.checkSat();
    if (r.isUnsat()) {
      throw std::runtime_error("Unsat");
    } else {
      for (const auto &relation : asserted) {
        Term term[2] = {relation_terms[relation][0],
                        relation_terms[relation][1]};
        auto lower = get_fpt64_value(solver.getValue(term[0]));
        auto upper = get_fpt64_value(solver.getValue(term[1]));
        relation_map[relation] = Probability(lower, upper);
      }
    }
  }
}
