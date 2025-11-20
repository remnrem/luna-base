#ifndef __LUNA_GAUSSIAN_HMM__
#define __LUNA_GAUSSIAN_HMM__

#include <stats/Eigen/Dense>
#include <vector>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <algorithm>

class gaussian_hmm_t {
public:
    gaussian_hmm_t(int n_states, int dim)
        : N_(n_states), M_(dim), debug_(false)
    {
        if (N_ <= 0 || M_ <= 0)
            throw std::runtime_error("N and M must be positive");

        pi_  = Eigen::VectorXd::Constant(N_, 1.0 / N_);
        A_   = Eigen::MatrixXd::Zero(N_, N_);
        mu_  = Eigen::MatrixXd::Zero(M_, N_);
        cov_.resize(N_);
        inv_cov_.resize(N_);
        log_det_cov_ = Eigen::VectorXd::Zero(N_);

        // Default: strong self-transitions, weak off-diagonals
        for (int i = 0; i < N_; ++i) {
            for (int j = 0; j < N_; ++j) {
                if (i == j) A_(i,j) = 0.90;
                else        A_(i,j) = 0.10 / double(N_ - 1);
            }
        }

        // Default emissions: zero-mean, identity covariance
        const double eps = 1e-8;
        for (int k = 0; k < N_; ++k) {
            cov_[k] = Eigen::MatrixXd::Identity(M_, M_);
            inv_cov_[k] = cov_[k].inverse();
            double det = cov_[k].determinant();
            if (det <= 0.0) det = eps;
            log_det_cov_(k) = std::log(det);
        }
    }

    int n_states() const { return N_; }
    int dim()      const { return M_; }

    void set_debug(bool debug) { debug_ = debug; }

    // Set initial state probabilities (will be normalized)
    void set_initial(const Eigen::VectorXd &pi) {
        if (pi.size() != N_)
            throw std::runtime_error("pi size mismatch");
        pi_ = pi;
        normalize_vec(pi_);
    }

    // Set transition matrix (rows will be normalized)
    void set_transition(const Eigen::MatrixXd &A) {
        if (A.rows() != N_ || A.cols() != N_)
            throw std::runtime_error("A size mismatch");
        A_ = A;
        for (int i = 0; i < N_; ++i)
            normalize_vec_row(A_, i);
    }

    // Set emission parameters: means (M x N) and full covariances (vector of M x M)
    void set_emission(const Eigen::MatrixXd &mu,
                      const std::vector<Eigen::MatrixXd> &cov)
    {
        if (mu.rows() != M_ || mu.cols() != N_)
            throw std::runtime_error("mu size mismatch");
        if ((int)cov.size() != N_)
            throw std::runtime_error("cov size mismatch");

        mu_ = mu;
        cov_ = cov;

        const double eps = 1e-8;
        for (int k = 0; k < N_; ++k) {
            if (cov_[k].rows() != M_ || cov_[k].cols() != M_)
                throw std::runtime_error("cov[k] size mismatch");

            // Regularize a bit, just in case
            cov_[k] += eps * Eigen::MatrixXd::Identity(M_, M_);
            inv_cov_[k] = cov_[k].inverse();
            double det = cov_[k].determinant();
            if (det <= 0.0) det = eps;
            log_det_cov_(k) = std::log(det);
        }
    }

    // Set emission parameters: means only; keep covariances as-is (or default)
    void set_emission(const Eigen::MatrixXd &mu)
    {
        if (mu.rows() != M_ || mu.cols() != N_)
            throw std::runtime_error("mu size mismatch");
        mu_ = mu;

        // If covariances not initialized, reset to identity
        if ((int)cov_.size() != N_) {
            cov_.assign(N_, Eigen::MatrixXd::Identity(M_, M_));
            inv_cov_.assign(N_, Eigen::MatrixXd::Identity(M_, M_));
            log_det_cov_.setZero();
        }
    }

    const Eigen::VectorXd&              pi()          const { return pi_; }
    const Eigen::MatrixXd&              A()           const { return A_; }
    const Eigen::MatrixXd&              mu()          const { return mu_; }
    const std::vector<Eigen::MatrixXd>& covariances() const { return cov_; }

    // ---------------- Training: single sequence ----------------

    // Train on a single contiguous sequence (no gaps).
    // obs[t] must be length M_.
    double train(const std::vector<Eigen::VectorXd> &obs,
                 int max_iters = 50,
                 double tol = 1e-4)
    {
        const int T = (int)obs.size();
        if (T < 2)
            throw std::runtime_error("Sequence too short for HMM training");
        for (int t = 0; t < T; ++t)
            if (obs[t].size() != M_)
                throw std::runtime_error("Observation dimension mismatch in train()");

        // Work arrays
        Eigen::MatrixXd alpha(T, N_);
        Eigen::MatrixXd beta (T, N_);
        Eigen::MatrixXd gamma(T, N_);
        std::vector<Eigen::MatrixXd> xi(T - 1, Eigen::MatrixXd(N_, N_));
        Eigen::VectorXd c(T);

        double prev_loglik = -std::numeric_limits<double>::infinity();

        for (int iter = 0; iter < max_iters; ++iter) {

            // E-step
            double loglik = forward_backward(obs, alpha, beta, gamma, xi, c);

            if (debug_) {
                std::cerr << "[train] iter " << iter
                          << " loglik = " << loglik << "\n";
            }

            // Non-finite guard
            if (!std::isfinite(loglik)) {
                if (debug_) {
                    std::cerr << "[train] non-finite loglik at iter "
                              << iter << ", breaking.\n";
                }
                break;
            }

            if (iter > 0) {
                double diff  = loglik - prev_loglik;
                double scale = std::max(std::abs(prev_loglik), 1.0);

                // If loglik decreased substantially, break
                if (diff < -1e-3 * scale) {
                    if (debug_) {
                        std::cerr << "[train] loglik decreased from "
                                  << prev_loglik << " to " << loglik
                                  << " at iter " << iter
                                  << ", breaking.\n";
                    }
                    break;
                }

                // Convergence: small positive improvement
                if (diff >= 0.0 && diff < scale * tol)
                    break;
            }
            prev_loglik = loglik;

            // M-step: accumulate sufficient stats from this single sequence
            Eigen::VectorXd gamma_sum     = Eigen::VectorXd::Zero(N_);
            Eigen::VectorXd gamma_sum_Tm1 = Eigen::VectorXd::Zero(N_);
            Eigen::MatrixXd xi_sum        = Eigen::MatrixXd::Zero(N_, N_);
            Eigen::MatrixXd num_mu        = Eigen::MatrixXd::Zero(M_, N_);
            std::vector<Eigen::MatrixXd> num_cov(N_, Eigen::MatrixXd::Zero(M_, M_));

            // Initial distribution
            pi_ = gamma.row(0).transpose();
            normalize_vec(pi_);

            // Gamma, xi, mu numerators
            for (int t = 0; t < T; ++t) {
                for (int i = 0; i < N_; ++i) {
                    double g = gamma(t, i);
                    gamma_sum(i) += g;
                    num_mu.col(i) += g * obs[t];
                    if (t < T - 1)
                        gamma_sum_Tm1(i) += g;
                }
                if (t < T - 1)
                    xi_sum += xi[t];
            }

            // Transition matrix
            for (int i = 0; i < N_; ++i) {
                double denom = gamma_sum_Tm1(i);
                const double eps = 1e-12;
                if (denom < eps) denom = eps;
                for (int j = 0; j < N_; ++j)
                    A_(i,j) = xi_sum(i,j) / denom;
                normalize_vec_row(A_, i);
            }

            // Means
            for (int i = 0; i < N_; ++i) {
                double denom = gamma_sum(i);
                const double eps = 1e-12;
                if (denom < eps) denom = eps;
                mu_.col(i) = num_mu.col(i) / denom;
            }

            // Covariances (use updated means)
            const double eps_cov = 1e-6;
            for (int t = 0; t < T; ++t) {
                for (int i = 0; i < N_; ++i) {
                    double g = gamma(t, i);
                    if (g <= 0.0) continue;
                    Eigen::VectorXd diff = obs[t] - mu_.col(i);
                    num_cov[i].noalias() += g * (diff * diff.transpose());
                }
            }

            for (int i = 0; i < N_; ++i) {
                double denom = gamma_sum(i);
                if (denom < eps_cov) denom = eps_cov;
                cov_[i] = num_cov[i] / denom;
                cov_[i].diagonal().array() += eps_cov;

                inv_cov_[i] = cov_[i].inverse();
                double det = cov_[i].determinant();
                if (det <= 0.0) det = eps_cov;
                log_det_cov_(i) = std::log(det);
            }
        }

        if (debug_) {
            std::cerr << "[train] final loglik = " << prev_loglik << "\n";
            debug_params();
        }

        return prev_loglik;
    }

    // ---------------- Training: multiple sequences ----------------

    // Train on multiple contiguous segments ("sequences"), e.g. recordings with gaps.
    // Each sequences[k] is an independent contiguous run.
    double train_multi(const std::vector<std::vector<Eigen::VectorXd>> &sequences,
                       int max_iters = 50,
                       double tol = 1e-4)
    {
        if (sequences.empty())
            throw std::runtime_error("No sequences provided to train_multi()");

        // Basic checks and total length
        std::size_t total_T = 0;
        for (const auto &seq : sequences) {
            if (seq.size() < 1) continue;
            for (const auto &x : seq)
                if (x.size() != M_)
                    throw std::runtime_error("Observation dimension mismatch in train_multi()");
            total_T += seq.size();
        }
        if (total_T < 2)
            throw std::runtime_error("Total length too short for HMM training");

        double prev_loglik = -std::numeric_limits<double>::infinity();

        for (int iter = 0; iter < max_iters; ++iter) {

            // Global accumulators across all sequences
            Eigen::VectorXd gamma_sum     = Eigen::VectorXd::Zero(N_);
            Eigen::VectorXd gamma_sum_Tm1 = Eigen::VectorXd::Zero(N_);
            Eigen::MatrixXd xi_sum        = Eigen::MatrixXd::Zero(N_, N_);
            Eigen::MatrixXd num_mu        = Eigen::MatrixXd::Zero(M_, N_);
            std::vector<Eigen::MatrixXd> num_cov(N_, Eigen::MatrixXd::Zero(M_, M_));
            Eigen::VectorXd pi_new        = Eigen::VectorXd::Zero(N_);

            double total_loglik = 0.0;

            // Store gammas per sequence for later covariance update
            std::vector<Eigen::MatrixXd> gamma_all(sequences.size());

            // E-step over each sequence
            for (std::size_t s = 0; s < sequences.size(); ++s) {

                const auto &obs = sequences[s];
                const int T = (int)obs.size();
                if (T < 2) {
                    gamma_all[s] = Eigen::MatrixXd(); // empty, skip
                    continue;
                }

                // Work arrays for this sequence
                Eigen::MatrixXd alpha(T, N_);
                Eigen::MatrixXd beta (T, N_);
                Eigen::MatrixXd gamma(T, N_);
                std::vector<Eigen::MatrixXd> xi(T - 1, Eigen::MatrixXd(N_, N_));
                Eigen::VectorXd c(T);

                double loglik = forward_backward(obs, alpha, beta, gamma, xi, c);
                total_loglik += loglik;

                gamma_all[s] = gamma;  // store for covariance update

                // Initial distribution contribution
                pi_new += gamma.row(0).transpose();

                // Accumulate gamma, xi, mu numerators
                for (int t = 0; t < T; ++t) {
                    for (int i = 0; i < N_; ++i) {
                        double g = gamma(t, i);
                        gamma_sum(i) += g;
                        num_mu.col(i) += g * obs[t];
                        if (t < T - 1)
                            gamma_sum_Tm1(i) += g;
                    }
                    if (t < T - 1)
                        xi_sum += xi[t];
                }
            }

            if (debug_) {
                std::cerr << "[train_multi] iter " << iter
                          << " total_loglik = " << total_loglik << "\n";
            }

            // Non-finite guard
            if (!std::isfinite(total_loglik)) {
                if (debug_) {
                    std::cerr << "[train_multi] non-finite total_loglik at iter "
                              << iter << ", breaking.\n";
                }
                break;
            }

            // Convergence check on total log-likelihood
            if (iter > 0) {
                double diff  = total_loglik - prev_loglik;
                double scale = std::max(std::abs(prev_loglik), 1.0);

                // If loglik decreased substantially, break
                if (diff < -1e-3 * scale) {
                    if (debug_) {
                        std::cerr << "[train_multi] total_loglik decreased from "
                                  << prev_loglik << " to " << total_loglik
                                  << " at iter " << iter
                                  << ", breaking.\n";
                    }
                    break;
                }

                // Convergence: small positive improvement
                if (diff >= 0.0 && diff < scale * tol)
                    break;
            }
            prev_loglik = total_loglik;

            // M-step using accumulated statistics

            // Initial distribution
            pi_ = pi_new;
            normalize_vec(pi_);

            // Transition matrix
            for (int i = 0; i < N_; ++i) {
                double denom = gamma_sum_Tm1(i);
                const double eps = 1e-12;
                if (denom < eps) denom = eps;
                for (int j = 0; j < N_; ++j)
                    A_(i,j) = xi_sum(i,j) / denom;
                normalize_vec_row(A_, i);
            }

            // Means
            for (int i = 0; i < N_; ++i) {
                double denom = gamma_sum(i);
                const double eps = 1e-12;
                if (denom < eps) denom = eps;
                mu_.col(i) = num_mu.col(i) / denom;
            }

            // Covariances (second pass, with updated means)
            const double eps_cov = 1e-2;
            for (std::size_t s = 0; s < sequences.size(); ++s) {
                const auto &obs = sequences[s];
                const int T = (int)obs.size();
                if (T < 2) continue;
                const auto &gamma = gamma_all[s];
                for (int t = 0; t < T; ++t) {
                    for (int i = 0; i < N_; ++i) {
                        double g = gamma(t, i);
                        if (g <= 0.0) continue;
                        Eigen::VectorXd diff = obs[t] - mu_.col(i);
                        num_cov[i].noalias() += g * (diff * diff.transpose());
                    }
                }
            }

            for (int i = 0; i < N_; ++i) {
                double denom = gamma_sum(i);
                if (denom < eps_cov) denom = eps_cov;
                cov_[i] = num_cov[i] / denom;
                // enforce floor on diagonal
                for (int d = 0; d < M_; ++d) {
                    if (cov_[i](d,d) < eps_cov)
                        cov_[i](d,d) = eps_cov;
                }

                inv_cov_[i] = cov_[i].inverse();
                double det = cov_[i].determinant();
                if (det <= 0.0) det = eps_cov;
                log_det_cov_(i) = std::log(det);
            }
        }

        if (debug_) {
            std::cerr << "[train_multi] final total_loglik = "
                      << prev_loglik << "\n";
            debug_params();
        }

        return prev_loglik;
    }

    // ---------------- Viterbi ----------------

    // Viterbi decoding for a single contiguous sequence.
    // Returns most likely state sequence (0..N-1).
    void viterbi(const std::vector<Eigen::VectorXd> &obs,
                 std::vector<int> &states) const
    {
        const int T = (int)obs.size();
        if (T == 0) {
            states.clear();
            return;
        }
        for (int t = 0; t < T; ++t)
            if (obs[t].size() != M_)
                throw std::runtime_error("Observation dimension mismatch in viterbi()");

        Eigen::MatrixXd delta(T, N_);
        Eigen::MatrixXi psi  (T, N_);

        // Initialization
        for (int i = 0; i < N_; ++i) {
            delta(0,i) = log_safe_prob(pi_(i)) +
                         gaussian_log_pdf_state(i, obs[0]);
            psi(0,i) = 0;
        }

        // Recursion
        for (int t = 1; t < T; ++t) {
            for (int j = 0; j < N_; ++j) {
                double best_val   = -std::numeric_limits<double>::infinity();
                int    best_state = 0;
                for (int i = 0; i < N_; ++i) {
                    double val = delta(t-1,i) + log_safe_prob(A_(i,j));
                    if (val > best_val) {
                        best_val   = val;
                        best_state = i;
                    }
                }
                delta(t,j) = best_val + gaussian_log_pdf_state(j, obs[t]);
                psi(t,j)   = best_state;
            }
        }

        // Backtrack
        states.assign(T, 0);
        int last_state = 0;
        double best = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < N_; ++i) {
            if (delta(T-1,i) > best) {
                best = delta(T-1,i);
                last_state = i;
            }
        }
        states[T-1] = last_state;
        for (int t = T-2; t >= 0; --t) {
            last_state = psi(t+1,last_state);
            states[t]  = last_state;
        }
    }

    // ---------------- Posteriors ----------------

    // Compute per-time posteriors gamma(t, i) = P(state i at t | observations).
    // gamma_out will be T x N_ (rows = time, cols = states).
    // Optionally returns log-likelihood if loglik != nullptr.
    void posteriors(const std::vector<Eigen::VectorXd> &obs,
                    Eigen::MatrixXd &gamma_out,
                    double *loglik = nullptr) const
    {
        const int T = (int)obs.size();
        if (T == 0)
            throw std::runtime_error("Empty sequence passed to posteriors()");
    
        for (int t = 0; t < T; ++t)
            if (obs[t].size() != M_)
                throw std::runtime_error("Observation dimension mismatch in posteriors()");
    
        Eigen::MatrixXd alpha(T, N_);
        Eigen::MatrixXd beta (T, N_);
        Eigen::MatrixXd gamma(T, N_);
        std::vector<Eigen::MatrixXd> xi;
        if (T > 1)
            xi.assign(T - 1, Eigen::MatrixXd(N_, N_));
        Eigen::VectorXd c(T);
    
        double ll = forward_backward(obs, alpha, beta, gamma, xi, c);
    
        gamma_out = gamma;

        if (debug_) {
            std::cerr << "[posteriors] loglik = " << ll << "\n";
            debug_gamma_stats(gamma);
            debug_state_occupancy(gamma);
            debug_gamma_rows(gamma, std::min(T, 10));
        }

        if (loglik) *loglik = ll;
    }

    // ---------------- Debug helpers (callable directly if needed) ----------------

    // Emission log-likelihoods B_log(t,i) = log p(obs_t | state i)
    void debug_emissions(const std::vector<Eigen::VectorXd> &obs,
                         Eigen::MatrixXd &B_log_out) const
    {
        const int T = (int)obs.size();
        B_log_out.resize(T, N_);
        for (int t = 0; t < T; ++t) {
            for (int i = 0; i < N_; ++i) {
                B_log_out(t,i) = gaussian_log_pdf_state(i, obs[t]);
            }
        }
    }

    // Entropy of gamma rows (how uniform the posteriors are)
    void debug_gamma_stats(const Eigen::MatrixXd &gamma) const
    {
        const int T = (int)gamma.rows();
        const int N = (int)gamma.cols();
        double min_H = std::numeric_limits<double>::infinity();
        double max_H = 0.0;

        for (int t = 0; t < T; ++t) {
            double H = 0.0;
            for (int i = 0; i < N; ++i) {
                double g = gamma(t,i);
                if (g > 0.0)
                    H -= g * std::log(g);
            }
            if (H < min_H) min_H = H;
            if (H > max_H) max_H = H;
        }
        std::cerr << "[debug_gamma_stats] entropy range: ["
                  << min_H << ", " << max_H
                  << "] vs log(N)=" << std::log((double)N_) << "\n";
    }

    // Per-state occupancy and argmax counts
    void debug_state_occupancy(const Eigen::MatrixXd &gamma) const
    {
        const int T = (int)gamma.rows();
        Eigen::VectorXd occ = gamma.colwise().sum();
        std::cerr << "[debug_state_occupancy] total gamma per state: "
                  << occ.transpose() << "\n";

        Eigen::VectorXi winners = Eigen::VectorXi::Zero(N_);
        for (int t = 0; t < T; ++t) {
            int argmax = 0;
            double best = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < N_; ++i) {
                double g = gamma(t,i);
                if (g > best) {
                    best = g;
                    argmax = i;
                }
            }
            winners(argmax)++;
        }
        std::cerr << "[debug_state_occupancy] argmax counts per state: "
                  << winners.transpose() << "\n";
    }

    // Print first T_print rows of gamma and their sums
    void debug_gamma_rows(const Eigen::MatrixXd &gamma,
                          int T_print = 10) const
    {
        int T = (int)gamma.rows();
        T_print = std::min(T_print, T);
        for (int t = 0; t < T_print; ++t) {
            double s = gamma.row(t).sum();
            std::cerr << "[debug_gamma_rows] t=" << t
                      << " sum=" << s
                      << " gamma=" << gamma.row(t) << "\n";
        }
    }

    // Print basic parameter sanity info
    void debug_params() const
    {
        std::cerr << "[debug_params] pi: " << pi_.transpose()
                  << " (sum=" << pi_.sum() << ")\n";

        std::cerr << "[debug_params] A row sums: ";
        for (int i = 0; i < N_; ++i) {
            double s = A_.row(i).sum();
            std::cerr << s << " ";
        }
        std::cerr << "\n";

        for (int i = 0; i < N_; ++i) {
            std::cerr << "[debug_params] state " << i
                      << " mu: " << mu_.col(i).transpose() << "\n";
            double min_diag = cov_[i].diagonal().minCoeff();
            double max_diag = cov_[i].diagonal().maxCoeff();
            std::cerr << "[debug_params] state " << i
                      << " cov diag range: "
                      << min_diag << " .. " << max_diag << "\n";
        }
    }

private:
    int N_; // number of states
    int M_; // observation dimension
    bool debug_;

    Eigen::VectorXd pi_;                   // N
    Eigen::MatrixXd A_;                    // N x N
    Eigen::MatrixXd mu_;                   // M x N (column i is mean of state i)
    std::vector<Eigen::MatrixXd> cov_;     // N of M x M
    std::vector<Eigen::MatrixXd> inv_cov_; // N of M x M
    Eigen::VectorXd log_det_cov_;          // N

    // log N(x | mu_i, cov_i)
    double gaussian_log_pdf_state(int k,
                                  const Eigen::VectorXd &x) const
    {
        const double two_pi = 6.28318530717958647692;
        Eigen::VectorXd diff = x - mu_.col(k);
        double quad = diff.transpose() * inv_cov_[k] * diff;
        double val = -0.5 * (M_ * std::log(two_pi) + log_det_cov_(k) + quad);

        if (!std::isfinite(val) && debug_) {
            std::cerr << "[gaussian_log_pdf_state] non-finite value for state "
                      << k << " quad=" << quad
                      << " log_det=" << log_det_cov_(k) << "\n";
        }
        return val;
    }

    static double log_safe_prob(double x) {
        if (x <= 0.0)
            return -std::numeric_limits<double>::infinity();
        return std::log(x);
    }
  
    static void normalize_vec(Eigen::VectorXd &v) {
        for (int i = 0; i < v.size(); ++i)
            if (!std::isfinite(v(i))) v(i) = 0.0;

        double s = v.sum();
        const double eps = 1e-12;
        if (s < eps) {
            v.setConstant(1.0 / v.size());
        } else {
            v /= s;
        }
    }

    static void normalize_vec_row(Eigen::MatrixXd &A, int i) {
        for (int j = 0; j < A.cols(); ++j)
            if (!std::isfinite(A(i,j))) A(i,j) = 0.0;

        double s = A.row(i).sum();
        const double eps = 1e-12;
        if (s < eps) {
            A.row(i).setConstant(1.0 / A.cols());
        } else {
            A.row(i) /= s;
        }
    }

    // Helper: log-sum-exp over a length-N vector
    static double log_sum_exp(const Eigen::VectorXd &v) {
        double maxv = v.maxCoeff();
        if (!std::isfinite(maxv)) return maxv; // all -inf etc.
        double sum = 0.0;
        for (int i = 0; i < v.size(); ++i)
            sum += std::exp(v(i) - maxv);
        return maxv + std::log(sum);
    }

    // Forward–backward entirely in log-space (no per-time scaling).
    //
    // - Returns log-likelihood log P(O | λ).
    // - Fills gamma(t,i) = P(q_t = i | O) in linear domain (each row sums to 1).
    // - Fills xi[t](i,j) = P(q_t = i, q_{t+1} = j | O) in linear domain.
    // - alpha, beta are returned in linear domain but are *not* used elsewhere.
    // - c(t) is set to 1.0 (only to keep the API).
    double forward_backward(
        const std::vector<Eigen::VectorXd> &obs,
        Eigen::MatrixXd &alpha,
        Eigen::MatrixXd &beta,
        Eigen::MatrixXd &gamma,
        std::vector<Eigen::MatrixXd> &xi,
        Eigen::VectorXd &c) const
    {
        const int T = (int)obs.size();
        if (T == 0)
            throw std::runtime_error("Empty sequence in forward_backward()");

        // Precompute log A
        Eigen::MatrixXd logA(N_, N_);
        for (int i = 0; i < N_; ++i)
            for (int j = 0; j < N_; ++j)
                logA(i,j) = log_safe_prob(A_(i,j));

        // Precompute log emissions: B_log(t,i) = log p(obs[t] | state i)
        Eigen::MatrixXd B_log(T, N_);
        for (int t = 0; t < T; ++t)
            for (int i = 0; i < N_; ++i)
                B_log(t,i) = gaussian_log_pdf_state(i, obs[t]);

        // Log-space forward and backward variables
        Eigen::MatrixXd log_alpha(T, N_);
        Eigen::MatrixXd log_beta (T, N_);

        // ---------- Forward ----------
        // t = 0
        for (int i = 0; i < N_; ++i)
            log_alpha(0,i) = log_safe_prob(pi_(i)) + B_log(0,i);

        // t = 1..T-1
        for (int t = 1; t < T; ++t) {
            for (int j = 0; j < N_; ++j) {
                Eigen::VectorXd tmp(N_);
                for (int i = 0; i < N_; ++i)
                    tmp(i) = log_alpha(t-1,i) + logA(i,j);
                double ls = log_sum_exp(tmp);
                log_alpha(t,j) = B_log(t,j) + ls;
            }
        }

        // Log-likelihood
        double loglik = log_sum_exp(log_alpha.row(T-1).transpose());
        if (!std::isfinite(loglik) && debug_) {
            std::cerr << "[forward_backward] non-finite loglik\n";
        }

        // ---------- Backward ----------
        // t = T-1
        for (int i = 0; i < N_; ++i)
            log_beta(T-1,i) = 0.0;  // log(1)

        // t = T-2..0
        for (int t = T-2; t >= 0; --t) {
            for (int i = 0; i < N_; ++i) {
                Eigen::VectorXd tmp(N_);
                for (int j = 0; j < N_; ++j)
                    tmp(j) = logA(i,j) + B_log(t+1,j) + log_beta(t+1,j);
                log_beta(t,i) = log_sum_exp(tmp);
            }
        }

        // ---------- Gamma ----------
        gamma.resize(T, N_);
        for (int t = 0; t < T; ++t) {
            Eigen::VectorXd log_gamma_t(N_);
            for (int i = 0; i < N_; ++i)
                log_gamma_t(i) = log_alpha(t,i) + log_beta(t,i) - loglik;

            // Normalize in linear domain
            double maxv = log_gamma_t.maxCoeff();
            Eigen::VectorXd tmp = (log_gamma_t.array() - maxv).exp();
            double s = tmp.sum();
            if (s <= 0.0) {
                // fallback: uniform if something pathological happens
                gamma.row(t).setConstant(1.0 / N_);
                if (debug_) {
                    std::cerr << "[forward_backward] gamma row " << t
                              << " had non-positive sum; set to uniform.\n";
                }
            } else {
                gamma.row(t) = (tmp / s).transpose();
            }
        }

        // ---------- Xi ----------
        if (T > 1 && !xi.empty()) {
            for (int t = 0; t < T - 1; ++t) {
                Eigen::MatrixXd &xitt = xi[t];
                xitt.resize(N_, N_);

                // compute log xi(i,j) ∝ log_alpha(t,i) + logA(i,j) + B_log(t+1,j) + log_beta(t+1,j) - loglik
                Eigen::VectorXd log_xi_vec(N_ * N_);
                int idx = 0;
                for (int i = 0; i < N_; ++i) {
                    for (int j = 0; j < N_; ++j) {
                        double val = log_alpha(t,i)
                                   + logA(i,j)
                                   + B_log(t+1,j)
                                   + log_beta(t+1,j)
                                   - loglik;
                        log_xi_vec(idx++) = val;
                    }
                }

                double maxv = log_xi_vec.maxCoeff();
                Eigen::VectorXd tmp = (log_xi_vec.array() - maxv).exp();
                double s = tmp.sum();
                if (s <= 0.0) {
                    xitt.setConstant(1.0 / double(N_ * N_));
                    if (debug_) {
                        std::cerr << "[forward_backward] xi at t=" << t
                                  << " had non-positive sum; set to uniform.\n";
                    }
                } else {
                    idx = 0;
                    for (int i = 0; i < N_; ++i)
                        for (int j = 0; j < N_; ++j)
                            xitt(i,j) = tmp(idx++) / s;
                }
            }
        }

        // ---------- Optionally return alpha, beta, c in linear domain ----------
        alpha.resize(T, N_);
        beta.resize(T, N_);
        c.resize(T);
        for (int t = 0; t < T; ++t) {
            // alpha: normalized for inspection/debugging
            double maxa = log_alpha.row(t).maxCoeff();
            Eigen::VectorXd tmpa = (log_alpha.row(t).array() - maxa).exp();
            double sa = tmpa.sum();
            if (sa <= 0.0) sa = 1.0;
            alpha.row(t) = (tmpa / sa).transpose();

            double maxb = log_beta.row(t).maxCoeff();
            Eigen::VectorXd tmpb = (log_beta.row(t).array() - maxb).exp();
            double sb = tmpb.sum();
            if (sb <= 0.0) sb = 1.0;
            beta.row(t) = (tmpb / sb).transpose();

            c(t) = 1.0;  // not used by callers
        }

        return loglik;
    }
};

#endif // __LUNA_GAUSSIAN_HMM__





