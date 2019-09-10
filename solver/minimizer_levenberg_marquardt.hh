#pragma once
#include "solver/base_solver.hh"
struct LevenbergMarguardtConfig{
  double min_norm_delta_x = 0.000000001;
  double tau = 0.000001;
  int max_continuous_fail_update = 10;
  double min_cost_function_gradient = 0.000000001;
};
namespace SLAMSolver {
class MinimizerLevenbergMarquardt {
public:
  MinimizerLevenbergMarquardt(std::shared_ptr<BaseSolver> solver_ptr, const LevenbergMarguardtConfig& config);
  bool solve(const int max_iterations);
private:
  double compute_initial_lambda();
  void add_lambda_to_hessian(const double lambda);
  void remove_lambda_from_hessian(const double lambda);
  double compute_gain_factor(const double cost_before_update, const double lambda);
private:
  std::shared_ptr<BaseSolver> solver_ptr_;
  LevenbergMarguardtConfig config_;
};
}