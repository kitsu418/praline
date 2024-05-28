#include "include/analysis.hpp"
#include <cstddef>
#include <iostream>
#include <optional>
#include <string>
#include <unistd.h>
#include <vector>

int main(int argc, char **argv) {
  int opt;
  bool is_legacy = false;
  std::optional<std::string> derivation_path, probability_path, output_path;

  auto print_helper_text = [&argv]() {
    std::cerr << "Usage: " << argv[0]
              << " -d DERIVATION_PATH -p PROBABILITY_PATH [-o OUTPUT_PATH]"
              << std::endl;
  };

  while ((opt = getopt(argc, argv, "hd:p:o:l")) != -1) {
    switch (opt) {
    case 'd':
      derivation_path = std::string(optarg);
      break;
    case 'p':
      probability_path = std::string(optarg);
      break;
    case 'o':
      output_path = std::string(optarg);
      break;
    case 'l':
      is_legacy = true;
      break;
    default:
      print_helper_text();
      return 1;
    }
  }

  if (!derivation_path.has_value() || !probability_path.has_value()) {
    print_helper_text();
    return 1;
  }

  auto analysis = new Analysis(derivation_path.value(),
                               probability_path.value(), is_legacy);
  if (is_legacy) {
    analysis->calculate_probability_legacy();
  } else {
    analysis->calculate_probability();
  }
  analysis->dump(output_path);
  delete analysis;

  return 0;
}
