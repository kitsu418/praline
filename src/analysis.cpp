#include "include/analysis.hpp"
#include "include/derivation.hpp"
#include "include/probability.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cvc5/cvc5.h>
#include <fstream>
#include <functional>
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

void Analysis::calculate_probability_legacy() {
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
            probability =
                probability.value().dep_disj(rule_probability_it->second);
          }
          continue;
        }

        Probability body_probability(1.0, 1.0);
        std::vector<size_t> unknown;

        for (size_t i = 0; i < body.size(); ++i) {
          auto relation_probability_it = relation_map.find(body[i]);
          if (relation_probability_it != relation_map.end()) {
            body_probability =
                body_probability.dep_conj(relation_probability_it->second);
          } else {
            unknown.push_back(i);
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
      }
      if (probability.has_value()) {
        relation_map[derivation.head] = probability.value();
      }
    }
  }
  solve_unknowns();
}

void Analysis::calculate_probability() {
  int node_num = 0;
  int edge_num = 0;
  for (const auto &derivation : derivations) {
    node_num++;
    for (auto &body : derivation.bodies) {
      edge_num += body.size();
    }
    // edge_num += derivation.bodies.size();
  }
  std::cerr << "node_num: " << node_num << std::endl;
  std::cerr << "edge_num: " << edge_num << std::endl;
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
            // probability = probability.value().dep_disj(body_probability);
            probability = operate(lhsd, probability.value(), lhsc,
                                  body_probability, kDisjunction);
            // std::cout << derivation.to_string() << std::endl;
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
  refine();
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
                   const std::string &probability_path, const bool &is_legacy) {
  this->derivations = Derivation::load(derivation_path);
  std::tie(this->relation_map, this->rule_map, this->dependent_class_map,
           this->dependent_map, this->queries) =
      Probability::load(probability_path, is_legacy);
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
    }

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

std::set<Relation> Analysis::dfs(Derivation &derivation,
                                 std::set<Derivation *> &visited,
                                 std::vector<Derivation *> &component) {
  if (visited.contains(&derivation)) {
    return {};
  }
  visited.insert(&derivation);
  std::set<Relation> leaves;

  if (relation_map.count(derivation.head)) {
    if (dependent_class_map.find(derivation.head) !=
        dependent_class_map.end()) {
      leaves.insert(derivation.head);
    }
  } else {
    std::optional<Probability> probability = Probability(0.0, 0.0);
    std::set<Relation> lhsd;

    for (const auto &body : derivation.bodies) {
      Probability body_probability(1.0, 1.0);
      std::set<Relation> lhsc;
      bool has_unknown = false;

      for (size_t i = 0; i < body.size(); ++i) {
        if (find(lhsc.begin(), lhsc.end(), body[i]) == lhsc.end()) {
          auto relation_probability = get_relation_probability(body[i]);
          if (!relation_probability.has_value()) {
            auto leaf = dfs(*derivations_index[body[i]], visited, component);
            leaves.merge(leaf);
          }
          if (relation_probability.has_value()) {
            for (int j = i + 1; j < body.size(); ++j) {
              auto condp_it =
                  dependent_map.find(std::make_pair(body[i], body[j]));
              if (condp_it != dependent_map.end()) {
                body_probability = body_probability.dep_conj(
                    condp_it->second.ind_conj(relation_probability.value()));
                lhsc.insert(body[j]);
                goto next_rel;
              }
            }
            body_probability =
                body_probability.dep_conj(relation_probability.value());
          next_rel:
            lhsc.insert(body[i]);
          } else {
            has_unknown = true;
            lhsc.insert(body[i]);
          }
        }
      }

      for (const auto &relation : body) {
        if (dependent_class_map.find(relation) != dependent_class_map.end()) {
          continue;
        }
        if (derivations_index.find(relation) != derivations_index.end()) {
          dfs(*derivations_index[relation], visited, component);
        }
      }
    }
  }

  component.push_back(&derivation);
  return leaves;
}

class Timer {
public:
  Timer() : start(std::chrono::high_resolution_clock::now()) {}
  ~Timer() {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cerr << "Time: " << duration.count() << "ms" << std::endl;
  }
  std::chrono::time_point<std::chrono::high_resolution_clock> start;
};

void Analysis::refine() {
  TermManager tm;
  Solver solver(tm);
  // solver.setOption("produce-models", "true");
  // solver.setOption("produce-unsat-cores", "true");
  solver.setOption("tlimit-per", "1000");
  solver.setLogic("ALL");

  Term rm = tm.mkRoundingMode(RoundingMode::ROUND_NEAREST_TIES_TO_EVEN);
  Sort fpt64 = tm.mkFloatingPointSort(11, 53);

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
  std::queue<Relation> roots;

  for (const auto &derivation : derivations) {
    if (queries.contains(derivation.head.name)) {
      roots.push(derivation.head);
    }
  }

  using state_t = uint64_t;
  const int DEPTH_LIMIT = 2;
  const double eps = 1e-8;

  std::function<void(const std::pair<Relation, int> &,
                     std::map<Relation, int> &, int &)>
      collect = [&](const std::pair<Relation, int> &current,
                    std::map<Relation, int> &collected, int &cnt) {
        auto relation = current.first;
        auto depth = current.second;

        if (depth == DEPTH_LIMIT ||
            (depth < DEPTH_LIMIT &&
             (!derivations_index.contains(relation) ||
              derivations_index[relation]->bodies.empty()))) {
          if (!collected.contains(relation)) {
            collected[relation] = cnt++;
          }
        } else if (depth < DEPTH_LIMIT) {
          for (const auto &body : derivations_index[relation]->bodies) {
            for (const auto &r : body) {
              collect({r, depth + 1}, collected, cnt);
            }
          }
        }
      };

  while (!roots.empty()) {
    auto root = roots.front();
    roots.pop();
    printf("%s\n", root.to_string().c_str());

    std::map<Relation, int> facts_id;
    int facts_num = 0;
    collect({root, 0}, facts_id, facts_num);
    solver.resetAssertions();

    // construct variables
    Term variables[(1 << facts_num)];
    for (state_t state = 0; state < (1ull << facts_num); ++state) {
      variables[state] = tm.mkConst(fpt64, "V_" + std::to_string(state));
      solver.assertFormula(
          tm.mkTerm(Kind::FLOATINGPOINT_LEQ, {variables[state], one}));
      solver.assertFormula(
          tm.mkTerm(Kind::FLOATINGPOINT_GEQ, {variables[state], zero}));
    }

    // add all sum to one constraint
    Term sum_to_one = variables[0];
    for (state_t state = 1; state < (1ull << facts_num); ++state) {
      sum_to_one = tm.mkTerm(Kind::FLOATINGPOINT_ADD,
                             {rm, sum_to_one, variables[state]});
    }
    solver.assertFormula(tm.mkTerm(Kind::FLOATINGPOINT_EQ, {sum_to_one, one}));

    std::function<std::map<state_t, double>(const std::pair<Relation, int> &)>
        add_constraints = [&](const std::pair<Relation, int> &current) {
          auto relation = current.first;
          auto depth = current.second;

          if (depth == DEPTH_LIMIT ||
              (depth < DEPTH_LIMIT &&
               (!derivations_index.contains(relation) ||
                derivations_index[relation]->bodies.empty()))) {
            // This is the fact case
            state_t mask = (1ull << facts_id[relation]);
            Term term = variables[mask];
            std::map<state_t, double> constraints{};
            constraints[mask] = 1.0;

            for (state_t state = mask + 1ull; state < (1ull << facts_num);
                 ++state) {
              if ((state & mask)) {
                term = tm.mkTerm(Kind::FLOATINGPOINT_ADD,
                                 {rm, term, variables[state]});
                constraints[state] = 1.0;
              }
            }

            solver.assertFormula(tm.mkTerm(
                Kind::FLOATINGPOINT_GEQ,
                {term, mk_fpt64_constant(relation_map[relation].lower_bound)}));
            solver.assertFormula(tm.mkTerm(
                Kind::FLOATINGPOINT_LEQ,
                {term, mk_fpt64_constant(relation_map[relation].upper_bound)}));

            return constraints;
          } else if (depth < DEPTH_LIMIT) {
            std::map<state_t, double> states;

            for (size_t i = 0; i < derivations_index[relation]->bodies.size();
                 ++i) {
              const auto &body = derivations_index[relation]->bodies[i];
              fprintf(stderr, "%s\n",
                      derivations_index[relation]->to_string().c_str());
              std::map<state_t, double> body_states;

              if (body.size() > 0) {
                body_states = add_constraints({body[0], depth + 1});
              }

              for (size_t j = 1; j < body.size(); ++j) {
                std::map<state_t, double> tmp1{};
                auto tmp2 = add_constraints({body[j], depth + 1});
                std::swap(tmp1, body_states);
                for (const auto &c1 : tmp1) {
                  for (const auto &c2 : tmp2) {
                    if (c1.first == c2.first) {
                      if (c1.second * c2.second < eps) {
                        body_states.erase(c1.first);
                      } else {
                        body_states[c1.first] = c1.second * c2.second;
                      }
                      break;
                    }
                  }
                }
              }

              auto rule_p = get_rule_probability(relation, body);
              if (rule_p.has_value()) {
                for (auto it = body_states.begin(); it != body_states.end();
                     ++it) {
                  it->second *= rule_p.value().lower_bound;
                }
              }

              for (const auto &c : body_states) {
                std::cerr << c.second << " * V_" << std::bitset<8>(c.first)
                          << " + ";
                // fprintf(stderr, "%.6lf * V_%lld + ", c.second, c.first);
              }
              // fprintf(stderr, "\n");
              std::cerr << std::endl;

              if (i == 0) {
                std::swap(states, body_states);
              } else {
                for (const auto &c : body_states) {
                  if (states.contains(c.first)) {
                    states[c.first] =
                        states[c.first] + c.second - states[c.first] * c.second;
                  } else {
                    states[c.first] = c.second;
                  }
                }
              }
            }

            // if (depth > 0) {
            //   auto it = states.begin();
            //   if (it != states.end()) {
            //     Term term = tm.mkTerm(
            //         Kind::FLOATINGPOINT_MULT,
            //         {rm, variables[it->first],
            //         mk_fpt64_constant(it->second)});
            //     for (++it; it != states.end(); ++it) {
            //       term =
            //           tm.mkTerm(Kind::FLOATINGPOINT_ADD,
            //                     {rm, term,
            //                      tm.mkTerm(Kind::FLOATINGPOINT_MULT,
            //                                {rm, variables[it->first],
            //                                 mk_fpt64_constant(it->second)})});
            //     }
            //     solver.assertFormula(tm.mkTerm(
            //         Kind::FLOATINGPOINT_GEQ,
            //         {term,
            //          mk_fpt64_constant(relation_map[relation].lower_bound)}));
            //     solver.assertFormula(tm.mkTerm(
            //         Kind::FLOATINGPOINT_LEQ,
            //         {term,
            //          mk_fpt64_constant(relation_map[relation].upper_bound)}));
            //   }
            // }

            return states;
          } else {
            return std::map<state_t, double>{};
          }
        };

    std::map<state_t, double> root_status;
    {
      Timer timer;
      root_status = add_constraints({root, 0});
    }
    std::cerr << root.to_string() << " = ";
    for (const auto &status : root_status) {
      std::cerr << status.second << " * V_" << std::bitset<8>(status.first)
                << " + ";
    }
    std::cerr << std::endl;

    if (root_status.empty()) {
      continue;
    }

    auto l = relation_map[root].lower_bound;
    auto u = relation_map[root].upper_bound;
    auto epsilon = (u - l) / 50.0;
    // auto epsilon = 0.00001;
    auto assertions = solver.getAssertions();
    // std::cerr << assertions << std::endl;

    auto it = root_status.begin();
    Term root_term =
        tm.mkTerm(Kind::FLOATINGPOINT_MULT,
                  {rm, variables[it->first], mk_fpt64_constant(it->second)});
    for (++it; it != root_status.end(); ++it) {
      root_term = tm.mkTerm(Kind::FLOATINGPOINT_ADD,
                            {rm, root_term,
                             tm.mkTerm(Kind::FLOATINGPOINT_MULT,
                                       {rm, variables[it->first],
                                        mk_fpt64_constant(it->second)})});
    }

    do {
      if (solver.getAssertions().empty()) {
        for (const auto &assertion : assertions) {
          solver.assertFormula(assertion);
        }
      }

      auto ll = l + epsilon > u ? u : l + epsilon;

      solver.assertFormula(tm.mkTerm(Kind::FLOATINGPOINT_LEQ,
                                     {root_term, mk_fpt64_constant(ll)}));
      solver.assertFormula(tm.mkTerm(Kind::FLOATINGPOINT_GEQ,
                                     {root_term, mk_fpt64_constant(l)}));

      if (solver.checkSat().isUnsat()) {
        l = ll;
        fprintf(stderr, "UNSAT: %.6lf\n", ll);
        solver.resetAssertions();
      } else {
        fprintf(stderr, "SAT: %.6lf\n", ll);
        break;
      }
    } while (l + epsilon < u);
    solver.resetAssertions();
    do {
      for (const auto &assertion : assertions) {
        solver.assertFormula(assertion);
      }

      auto uu = u - epsilon < l ? l : u - epsilon;

      solver.assertFormula(tm.mkTerm(Kind::FLOATINGPOINT_LEQ,
                                     {root_term, mk_fpt64_constant(u)}));
      solver.assertFormula(tm.mkTerm(Kind::FLOATINGPOINT_GEQ,
                                     {root_term, mk_fpt64_constant(uu)}));

      if (solver.checkSat().isUnsat()) {
        u = uu;
        fprintf(stderr, "UNSAT: %.6lf\n", uu);
        solver.resetAssertions();
      } else {
        fprintf(stderr, "SAT: %.6lf\n", uu);
        break;
      }
    } while (u - epsilon > l);
    relation_map[root] = Probability(l, u);
  }
}