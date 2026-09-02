//    --------------------------------------------------------------------
//    Statistical evaluation for subject-level outer-OOF predictions.
//    --------------------------------------------------------------------

//    This file is part of Luna.

//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    Luna is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with Luna. If not, see <http://www.gnu.org/licenses/>.

//    Please see LICENSE.txt for more details.

//    --------------------------------------------------------------------

#include "stats/dpp-evaluation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>

namespace {
const double NaN = std::numeric_limits<double>::quiet_NaN();

void check_same_size(const std::vector<double> &a,
                     const std::vector<double> &b) {
  if (a.size() != b.size())
    throw std::invalid_argument(
        "paired evaluation vectors must have equal length");
}

std::vector<std::pair<double, double>>
finite_pairs(const std::vector<double> &observed,
             const std::vector<double> &predicted) {
  check_same_size(observed, predicted);
  std::vector<std::pair<double, double>> result;
  for (std::size_t i = 0; i < observed.size(); ++i)
    if (std::isfinite(observed[i]) && std::isfinite(predicted[i]))
      result.push_back({observed[i], predicted[i]});
  return result;
}

std::vector<double> finite_differences(const std::vector<double> &a,
                                       const std::vector<double> &b) {
  check_same_size(a, b);
  std::vector<double> differences;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (std::isfinite(a[i]) && std::isfinite(b[i]))
      differences.push_back(a[i] - b[i]);
  return differences;
}

double regularized_beta(double a, double b, double x) {
  if (x <= 0.0)
    return 0.0;
  if (x >= 1.0)
    return 1.0;
  const int max_iterations = 200;
  const double epsilon = 3.0e-14;
  const double tiny = 1.0e-300;
  const double log_factor = std::lgamma(a + b) - std::lgamma(a) -
                            std::lgamma(b) + a * std::log(x) +
                            b * std::log(1.0 - x);
  const double factor = std::exp(log_factor);
  double qab = a + b;
  double qap = a + 1.0;
  double qam = a - 1.0;
  double c = 1.0;
  double d = 1.0 - qab * x / qap;
  d = std::fabs(d) < tiny ? tiny : d;
  d = 1.0 / d;
  double h = d;
  for (int m = 1; m <= max_iterations; ++m) {
    const double m2 = 2.0 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    d = std::fabs(d) < tiny ? tiny : d;
    c = 1.0 + aa / c;
    c = std::fabs(c) < tiny ? tiny : c;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    d = std::fabs(d) < tiny ? tiny : d;
    c = 1.0 + aa / c;
    c = std::fabs(c) < tiny ? tiny : c;
    d = 1.0 / d;
    const double delta = d * c;
    h *= delta;
    if (std::fabs(delta - 1.0) < epsilon)
      break;
  }
  return factor * h / a;
}

double student_t_two_sided_p(double statistic, std::size_t n) {
  if (n < 2)
    return 1.0;
  if (std::isinf(statistic))
    return 0.0;
  if (!std::isfinite(statistic))
    return 1.0;
  const double degrees_of_freedom = static_cast<double>(n - 1);
  const double t = std::fabs(statistic);
  const double x = degrees_of_freedom / (degrees_of_freedom + t * t);
  const double tail = 0.5 * regularized_beta(degrees_of_freedom / 2.0, 0.5, x);
  return std::max(0.0, std::min(1.0, 2.0 * tail));
}
} // namespace

namespace dpp_evaluation {
double rmse(const std::vector<double> &observed,
            const std::vector<double> &predicted) {
  const auto pairs = finite_pairs(observed, predicted);
  if (pairs.empty())
    return NaN;
  double squared_error = 0.0;
  for (const auto &pair : pairs)
    squared_error += (pair.first - pair.second) * (pair.first - pair.second);
  return std::sqrt(squared_error / pairs.size());
}

double r2(const std::vector<double> &observed,
          const std::vector<double> &predicted) {
  const auto pairs = finite_pairs(observed, predicted);
  if (pairs.empty())
    return NaN;
  double mean = 0.0;
  for (const auto &pair : pairs)
    mean += pair.first;
  mean /= pairs.size();
  double total = 0.0;
  double residual = 0.0;
  for (const auto &pair : pairs) {
    total += (pair.first - mean) * (pair.first - mean);
    residual += (pair.first - pair.second) * (pair.first - pair.second);
  }
  return total > 0.0 ? 1.0 - residual / total : NaN;
}

double brier_score(const std::vector<double> &observed,
                   const std::vector<double> &probabilities) {
  const auto pairs = finite_pairs(observed, probabilities);
  if (pairs.empty())
    return NaN;
  double total = 0.0;
  for (const auto &pair : pairs)
    total += (pair.first - pair.second) * (pair.first - pair.second);
  return total / pairs.size();
}

double log_loss(const std::vector<double> &observed,
                const std::vector<double> &probabilities) {
  const auto pairs = finite_pairs(observed, probabilities);
  if (pairs.empty())
    return NaN;
  double total = 0.0;
  for (const auto &pair : pairs) {
    const double probability =
        std::max(1.0e-15, std::min(1.0 - 1.0e-15, pair.second));
    total -= pair.first * std::log(probability) +
             (1.0 - pair.first) * std::log(1.0 - probability);
  }
  return total / pairs.size();
}

double auc(const std::vector<double> &observed,
           const std::vector<double> &scores) {
  const auto pairs = finite_pairs(observed, scores);
  std::vector<double> positive_scores;
  std::vector<double> negative_scores;
  for (const auto &pair : pairs) {
    if (pair.first == 1.0)
      positive_scores.push_back(pair.second);
    else if (pair.first == 0.0)
      negative_scores.push_back(pair.second);
    else
      throw std::invalid_argument("AUC outcomes must be coded 0/1");
  }
  if (positive_scores.empty() || negative_scores.empty())
    return NaN;
  double wins = 0.0;
  for (double positive : positive_scores)
    for (double negative : negative_scores)
      if (positive > negative)
        wins += 1.0;
      else if (positive == negative)
        wins += 0.5;
  return wins / (positive_scores.size() * negative_scores.size());
}

prediction_metrics_t evaluate_predictions(const std::vector<double> &observed,
                                          const std::vector<double> &predicted,
                                          bool binary_outcome) {
  const auto pairs = finite_pairs(observed, predicted);
  prediction_metrics_t result;
  result.rmse = result.r2 = result.brier = result.log_loss = result.auc = NaN;
  result.n = pairs.size();
  if (result.n == 0)
    return result;
  result.rmse = rmse(observed, predicted);
  if (binary_outcome) {
    result.brier = brier_score(observed, predicted);
    result.log_loss = log_loss(observed, predicted);
    result.auc = auc(observed, predicted);
  } else
    result.r2 = r2(observed, predicted);
  return result;
}

double paired_mean_loss_difference(const std::vector<double> &loss_a,
                                   const std::vector<double> &loss_b) {
  const auto differences = finite_differences(loss_a, loss_b);
  if (differences.empty())
    return NaN;
  double total = 0.0;
  for (double difference : differences)
    total += difference;
  return total / differences.size();
}

paired_test_result_t paired_t_test(const std::vector<double> &loss_a,
                                   const std::vector<double> &loss_b) {
  const auto differences = finite_differences(loss_a, loss_b);
  paired_test_result_t result;
  result.n = differences.size();
  if (differences.empty()) {
    result.mean_difference = result.standard_deviation = NaN;
    return result;
  }
  for (double difference : differences)
    result.mean_difference += difference;
  result.mean_difference /= result.n;
  if (result.n < 2) {
    result.standard_deviation = NaN;
    result.statistic = NaN;
    return result;
  }
  double sum_squared = 0.0;
  for (double difference : differences)
    sum_squared += (difference - result.mean_difference) *
                   (difference - result.mean_difference);
  result.standard_deviation = std::sqrt(sum_squared / (result.n - 1));
  if (result.standard_deviation == 0.0)
    result.statistic = result.mean_difference == 0.0
                           ? 0.0
                           : (result.mean_difference > 0.0
                                  ? std::numeric_limits<double>::infinity()
                                  : -std::numeric_limits<double>::infinity());
  else
    result.statistic =
        result.mean_difference /
        (result.standard_deviation / std::sqrt(static_cast<double>(result.n)));
  result.p_value = student_t_two_sided_p(result.statistic, result.n);
  return result;
}

permutation_test_result_t
paired_permutation_test(const std::vector<double> &loss_a,
                        const std::vector<double> &loss_b,
                        std::size_t permutations, std::uint64_t seed) {
  const auto differences = finite_differences(loss_a, loss_b);
  permutation_test_result_t result;
  result.n = differences.size();
  result.permutations = permutations;
  if (differences.empty() || permutations == 0) {
    result.observed_difference = differences.empty() ? NaN : 0.0;
    result.p_value = differences.empty() ? NaN : 1.0;
    return result;
  }
  result.observed_difference = 0.0;
  for (double difference : differences)
    result.observed_difference += difference;
  result.observed_difference /= differences.size();
  std::mt19937_64 generator(seed);
  std::bernoulli_distribution sign(0.5);
  const double threshold = std::fabs(result.observed_difference);
  for (std::size_t permutation = 0; permutation < permutations; ++permutation) {
    double permuted = 0.0;
    for (double difference : differences)
      permuted += sign(generator) ? difference : -difference;
    permuted /= differences.size();
    if (std::fabs(permuted) >= threshold)
      ++result.extreme_permutations;
  }
  result.p_value = static_cast<double>(result.extreme_permutations + 1) /
                   static_cast<double>(permutations + 1);
  return result;
}

auc_chance_test_result_t
auc_chance_permutation_test(const std::vector<double> &observed,
                            const std::vector<double> &scores,
                            std::size_t permutations, std::uint64_t seed) {
  const auto pairs = finite_pairs(observed, scores);
  auc_chance_test_result_t result;
  result.n = pairs.size();
  result.permutations = permutations;
  if (pairs.empty() || permutations == 0) {
    result.auc = pairs.empty() ? NaN : auc(observed, scores);
    result.difference_from_chance = pairs.empty() ? NaN : result.auc - 0.5;
    result.p_value = pairs.empty() ? NaN : 1.0;
    return result;
  }

  std::vector<double> finite_observed;
  std::vector<double> finite_scores;
  for (const auto &pair : pairs) {
    finite_observed.push_back(pair.first);
    finite_scores.push_back(pair.second);
  }
  result.auc = auc(finite_observed, finite_scores);
  result.difference_from_chance = result.auc - 0.5;
  std::mt19937_64 generator(seed);
  for (std::size_t permutation = 0; permutation < permutations; ++permutation) {
    std::shuffle(finite_scores.begin(), finite_scores.end(), generator);
    const double permuted_auc = auc(finite_observed, finite_scores);
    if (std::fabs(permuted_auc - 0.5) >=
        std::fabs(result.difference_from_chance))
      ++result.extreme_permutations;
  }
  result.p_value = static_cast<double>(result.extreme_permutations + 1) /
                   static_cast<double>(permutations + 1);
  return result;
}
} // namespace dpp_evaluation
