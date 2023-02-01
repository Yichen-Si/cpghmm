#ifndef CTHMM_H
#define CTHMM_H

#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <cmath>
#include <cfloat>
#include <iostream>

#include "brent_obj.h"

#define OPTIM_ENABLE_EIGEN_WRAPPERS
// #include "optim/optim.hpp"
#include "optim.hpp"

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXXd;
using Eigen::ArrayXd;

struct History {
    int32_t iter;
    int32_t failed;
    double loo_post, loo_map, viterbi_ll;
    ArrayXd theta;
};

struct Box {
    VectorXd upper, lower;
};

class cthmm {
public:
    int32_t n_state, n_obs, n_category, record_freq;
    std::vector<int32_t>& obs;
    std::vector<float>& distance;
    double dist_scale, max_scale, min_scale;
    ArrayXXd Amtx, Emtx; // Transition (A) and Emissition (E) probabilities
    ArrayXXd alpha, beta, loo, marginal;
    ArrayXd theta, init_state_prob, max_ll;
    Eigen::Array<int32_t,Eigen::Dynamic,Eigen::Dynamic> phi;
    std::vector<int32_t> viterbi_path;
    double loo_post, loo_map, viterbi_ll;
    ArrayXd update_size = ArrayXd::Zero(3);
    std::vector<History*> track;

    cthmm(std::vector<int32_t>& _obs, std::vector<float>& _dist, double _s,
          std::vector<double>& _scale, ArrayXXd _A, ArrayXXd _E, ArrayXd _pi, int32_t _rec = 5) :
          obs(_obs), distance(_dist), dist_scale(_s), Amtx(_A), Emtx(_E), init_state_prob(_pi)
          {
              record_freq= _rec;
              n_state    = _scale.size();
              n_obs      = obs.size();
              n_category = Emtx.cols();
              alpha = MatrixXd::Zero(n_state, n_obs);
              beta  = MatrixXd::Zero(n_state, n_obs);
              theta.resize(n_state);
              for (int32_t i = 0; i < n_state; ++i) {
                  theta(i) = 1./_scale[i];
              }
              min_scale = 20 * dist_scale;
              max_scale = 3e6* dist_scale;
          }
    ~cthmm() {
        for (auto & v : track) {delete v;}
    }

    // Basic HMM functions
    void forward();
    void backward();
    void conditional_prob();
    void leave_one_out();
    void leave_one_out_composite_posterior();
    void leave_one_out_composite_map();
    void viterbi();
    // Compute weighted average for emission and transition probabilities
    void update_matrix();

    // Optimize LOO MAP likelihood. Update individual transition rate iteratively
    void mixed_optim( int32_t max_iter = 20, int32_t max_iter_inner = 20, double tol = 1e-8,  double tol_inner = 1e-3, int32_t criterion = 0);

    // EM with Brent's method for transition rate
    void EM( int32_t max_iter_EM = 20, int32_t max_iter_inner = 20, double tol_EM = 1e-8, double tol_inner = 1e-8, int32_t optim_method = 0);
    void min_obj_theta(ArrayXd x, ArrayXd& res);
    void optim_brent_theta(int32_t max_iter, double tol, ArrayXd& arg);
    double min_obj_theta_indivisual(int32_t state_idx, double x);
    double optim_brent_theta_individual(int32_t state_idx, int32_t max_iter, double tol);

    // This is not used anymore, but derivative computation can be found here
    void NewtonRaphson(int32_t max_iter_inner, double tol_inner, ArrayXd& theta_new);
};

// Functions for using optim library to update \theta
double obj_map_wrap(const VectorXd& vals_inp, VectorXd* grad_out, void* opt_data);
void nm_loo_map_optim(cthmm* hmm_obj, int32_t max_iter, int32_t max_iter_inner, double tol, double bound_scale_lower=0.01, double bound_scale_upper=5e3);
VectorXd constr_fn(const VectorXd& vals_inp, MatrixXd* jacob_out, void* constr_data);

// Function for EM update with fixed transition rate
void em_fix_transition_rate(cthmm* hmm_obj, int32_t max_iter, double tol);

#endif
