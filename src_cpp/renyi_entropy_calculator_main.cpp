
#include <iostream>
#include <ranges>

#include <boost/program_options.hpp>

#include "random_samples.h"
#include "renyi_entropy.h"
#include "utils.h"

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  std::string directory;
  unsigned int maximal_neighborhood;
  std::vector<std::vector<unsigned int>> history_firsts, future_firsts,
      history_seconds;
  std::vector<double> alpha;
  std::string alphas_string, history_first_string, future_first_string,
      history_second_string;
  po::options_description desc("Allowed options");
  desc.add_options()("help, h", "produce help message")(
      "directory, d", po::value<std::string>(), "Folder to export results")(
      "history_first", po::value<std::string>()->composing(),
      "History of the first timeries")("future_first",
                                       po::value<std::string>()->composing(),
                                       "Future of the first timeries")(
      "history_second", po::value<std::string>()->composing(),
      "History of the second timeseries")(
      "maximal_neighborhood", po::value<unsigned int>(),
      "Maximal neighborhood")("alpha", po::value<std::string>()->composing(),
                              "Renyi entropy parameter");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  if (vm.count("directory")) {
    directory = vm["directory"].as<std::string>();
  }
  if (vm.count("history_first")) {
    history_first_string = vm["history_first"].as<std::string>();
    chop_string_arrays<std::string, unsigned int>(
        history_first_string, ',', history_firsts,
        [](const std::string str, std::size_t *pos) -> unsigned int {
          return std::stoul(str, pos, 10);
        });
  } else {
    std::cerr << "history_first option must be used" << std::endl;
  }
  if (vm.count("future_first")) {
    future_first_string = vm["future_first"].as<std::string>();
    chop_string_arrays<std::string, unsigned int>(
        future_first_string, ',', future_firsts,
        [](const std::string str, std::size_t *pos) -> unsigned int {
          return std::stoul(str, pos, 10);
        });
  } else {
    std::cerr << "future_first option must be used" << std::endl;
  }
  if (vm.count("history_second")) {
    history_second_string = vm["history_second"].as<std::string>();
    chop_string_arrays<std::string, unsigned int>(
        history_second_string, ',', history_seconds,
        [](const std::string str, std::size_t *pos) -> unsigned int {
          return std::stoul(str, pos, 10);
        });
  } else {
    std::cerr << "history_second option must be used" << std::endl;
  }
  if (vm.count("maximal_neighborhood")) {
    maximal_neighborhood = vm["maximal_neighborhood"].as<unsigned int>();
  }
  if (vm.count("alpha")) {
    alphas_string = vm["alpha"].as<std::string>();
    convert_array<std::string, double>(alphas_string, alpha, std::stod);
  }
  unsigned int postselection_X_future;
  unsigned int postselection_X_history;
  unsigned int postselection_Y_history;

  for (auto &item : alpha) {
    std::cout << item << std::endl;
  }
  Eigen::MatrixXd dataset;
  random_samples::samples_normal_distribution_uncorrelated(dataset, 0, 1, 100,
                                                           2);

  for (auto swap_datasets : {true, false}) {
    for (auto shuffle_dataset : {true, false}) {
      auto [dataset1, dataset2] =
          renyi_entropy::renyi_entropy<double>::prepare_dataset(
              dataset, swap_datasets, shuffle_dataset, 1, 1);

      for (auto &future_first : future_firsts) {
        for (auto &histories_first : history_firsts) {
          for (auto &histories_second : history_seconds) {
            std::map<std::string, std::any> configuration_prepare_dataset = {
                std::pair<std::string, std::any>("transpose", true),
                std::pair<std::string, std::any>("history_index_x",
                                                 histories_first),
                std::pair<std::string, std::any>("history_index_y",
                                                 histories_second),
                std::pair<std::string, std::any>("future_index_x",
                                                 future_first),
                std::pair<std::string, std::any>("postselection_y_fut",
                                                 postselection_X_future),
                std::pair<std::string, std::any>("postselection_z_hist",
                                                 postselection_Y_history),
                std::pair<std::string, std::any>("postselection_y_hist",
                                                 postselection_X_history)};
            auto [y_future, y_history, z_history] = renyi_entropy::
                renyi_entropy<double>::PrepareDatasetForTransferEntropy(
                    dataset1, dataset2, configuration_prepare_dataset);

            auto indices_to_use =
                std::ranges::iota_view{1U, maximal_neighborhood};
            std::map<std::string, std::any> configuration_renyi_entropy = {
                std::pair<std::string, std::any>("transpose", true),
                std::pair<std::string, std::any>("axis_to_join", 0),
                std::pair<std::string, std::any>("method", "LeonenkoProzanto"),
                std::pair<std::string, std::any>("alphas", alpha),
                std::pair<std::string, std::any>("enhanced_calculation", true),
                std::pair<std::string, std::any>("indices_to_use",
                                                 indices_to_use)};

            auto transfer_entropy = renyi_entropy::renyi_entropy<double>::
                renyi_conditional_information_transfer(
                    y_future, y_history, z_history,
                    configuration_renyi_entropy);
          }
        }
      }
    }
  }
}
