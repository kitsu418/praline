#include "include/analysis.hpp"
#include <cstddef>
#include <iostream>
#include <optional>
#include <string>
#include <unistd.h>
#include <vector>

int main(int argc, char **argv) {
  int opt;
  std::optional<std::string> derivation_path, probability_path, output_path;

  auto print_helper_text = [&argv]() {
    std::cerr << "Usage: " << argv[0]
              << " -d DERIVATION_PATH -p PROBABILITY_PATH [-o OUTPUT_PATH]"
              << std::endl;
  };

  while ((opt = getopt(argc, argv, "hd:p:o:")) != -1) {
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
    default:
      print_helper_text();
      return 1;
    }
  }

  if (!derivation_path.has_value() || !probability_path.has_value()) {
    print_helper_text();
    return 1;
  }

  auto analysis =
      new Analysis(derivation_path.value(), probability_path.value());
  analysis->calculate_probability();
  analysis->dump(output_path);
  delete analysis;

  return 0;
}
