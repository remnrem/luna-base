
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

#ifndef __KMEANS_EIGEN_H__
#define __KMEANS_EIGEN_H__

#include <stats/Eigen/Dense>
#include <vector>
#include <random>
#include <limits>
#include <algorithm>

struct kmeans_result_t {
    Eigen::MatrixXd centroids;  // K x D
    Eigen::VectorXi labels;     // N
};

// X: N x D (rows = samples)
// K: number of clusters
// max_iters: max Lloyd iterations
// tol: convergence tolerance on centroid movement (L2)
// seed: RNG seed for reproducible init
inline kmeans_result_t kmeans(const Eigen::MatrixXd &X,
                              int K,
                              int max_iters = 100,
                              double tol = 1e-4,
                              unsigned int seed = 12345)
{
    const int N = static_cast<int>(X.rows());
    const int D = static_cast<int>(X.cols());
    if (N == 0 || D == 0 || K <= 0)
        throw std::runtime_error("kmeans: invalid dimensions");
    if (K > N)
        throw std::runtime_error("kmeans: K > N not supported");

    kmeans_result_t res;
    res.centroids.resize(K, D);
    res.labels.setZero(N);

    // ----- 1) Randomly pick K distinct points as initial centroids -----
    std::mt19937 rng(seed);
    std::vector<int> indices(N);
    for (int i = 0; i < N; ++i) indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), rng);
    for (int k = 0; k < K; ++k)
        res.centroids.row(k) = X.row(indices[k]);

    Eigen::MatrixXd old_centroids = res.centroids;

    // ----- 2) Lloyd iterations -----
    for (int iter = 0; iter < max_iters; ++iter) {

        // Assignment step
        for (int n = 0; n < N; ++n) {
            double best_dist = std::numeric_limits<double>::infinity();
            int best_k = 0;
            for (int k = 0; k < K; ++k) {
                double dist = (X.row(n) - res.centroids.row(k)).squaredNorm();
                if (dist < best_dist) {
                    best_dist = dist;
                    best_k = k;
                }
            }
            res.labels(n) = best_k;
        }

        // Update step
        res.centroids.setZero();
        Eigen::VectorXi counts = Eigen::VectorXi::Zero(K);
        for (int n = 0; n < N; ++n) {
            int k = res.labels(n);
            res.centroids.row(k) += X.row(n);
            counts(k)++;
        }

        // Handle empty clusters: reinit to random point
        for (int k = 0; k < K; ++k) {
            if (counts(k) > 0) {
                res.centroids.row(k) /= static_cast<double>(counts(k));
            } else {
                // Reinitialize to a random data point
                std::uniform_int_distribution<int> uni(0, N - 1);
                int idx = uni(rng);
                res.centroids.row(k) = X.row(idx);
            }
        }

        // Check convergence (movement of centroids)
        double move = (res.centroids - old_centroids).norm();
        if (move < tol)
            break;

        old_centroids = res.centroids;
    }

    return res;
}


#endif 
