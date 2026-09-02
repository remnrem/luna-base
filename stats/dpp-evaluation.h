//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Luna is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Luna. If not, see <http://www.gnu.org/licenses/>.
//
//    Please see LICENSE.txt for more details.
//
//    --------------------------------------------------------------------
//
//    Statistical evaluation for subject-level outer-OOF predictions.
//    --------------------------------------------------------------------

#ifndef __LUNA_DPP_EVALUATION_H__
#define __LUNA_DPP_EVALUATION_H__

#include <cstddef>
#include <cstdint>
#include <vector>

namespace dpp_evaluation {
struct prediction_metrics_t {
  std::size_t n = 0;
  double rmse = 0.0;
  double r2 = 0.0;
  double brier = 0.0;
  double log_loss = 0.0;
  double auc = 0.0;
};

struct paired_test_result_t {
  std::size_t n = 0;
  double mean_difference = 0.0;
  double standard_deviation = 0.0;
  double statistic = 0.0;
  double p_value = 1.0;
};

struct permutation_test_result_t {
  std::size_t n = 0;
  std::size_t permutations = 0;
  std::size_t extreme_permutations = 0;
  double observed_difference = 0.0;
  double p_value = 1.0;
};

struct auc_chance_test_result_t {
  std::size_t n = 0;
  double auc = 0.0;
  double difference_from_chance = 0.0;
  std::size_t permutations = 0;
  std::size_t extreme_permutations = 0;
  double p_value = 1.0;
};

double rmse(const std::vector<double> &observed,
            const std::vector<double> &predicted);
double r2(const std::vector<double> &observed,
          const std::vector<double> &predicted);
double brier_score(const std::vector<double> &observed,
                   const std::vector<double> &probabilities);
double log_loss(const std::vector<double> &observed,
                const std::vector<double> &probabilities);
double auc(const std::vector<double> &observed,
           const std::vector<double> &scores);

prediction_metrics_t evaluate_predictions(const std::vector<double> &observed,
                                          const std::vector<double> &predicted,
                                          bool binary_outcome);

double paired_mean_loss_difference(const std::vector<double> &loss_a,
                                   const std::vector<double> &loss_b);
paired_test_result_t paired_t_test(const std::vector<double> &loss_a,
                                   const std::vector<double> &loss_b);
permutation_test_result_t paired_permutation_test(
    const std::vector<double> &loss_a, const std::vector<double> &loss_b,
    std::size_t permutations = 10000, std::uint64_t seed = 1);

auc_chance_test_result_t auc_chance_permutation_test(
    const std::vector<double> &observed, const std::vector<double> &scores,
    std::size_t permutations = 10000, std::uint64_t seed = 1);
} // namespace dpp_evaluation

#endif
