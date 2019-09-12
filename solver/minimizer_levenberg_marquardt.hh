#pragma once
#include "solver/base_solver.hh"
struct LevenbergMarguardtConfig{
  double min_norm_delta_x = 0.000000001;
  double tau = 0.000001;
  int max_continuous_fail_update = 10;
  double min_cost_function_gradient = 0.000000001;
};
namespace SLAMSolver {
/*!
 * @brief Non-linear least square minimizer with Levenberg-Marguardt method
 */
class MinimizerLevenbergMarquardt {
public:
  /*!
   * @brief Constructor
   * @param solver_ptr A solver pointer containing the problem
   * @param config Configuration of minimizer
   */
  MinimizerLevenbergMarquardt(std::shared_ptr<BaseSolver> solver_ptr, const LevenbergMarguardtConfig& config);
  /*!
   * @brief Minimize the non-linear least square cost
   * @param max_iterations maximum number of iterations
   * @return
   */
  bool minimize(const int max_iterations);
private:
  /*!
   * @brief Compute initial lambda value used in LM
   * @return lambda value
   */
  double compute_initial_lambda();
  /*!
   * @brief Add lambda to hessian matrix in solver_ptr
   * @param lambda lambda value hessian += lambda * Identity
   */
  void add_lambda_to_hessian(const double lambda);
  /*!
   * @brief Remove added lambda from hessian matrix in solver_ptr
   * @param lambda Lambda value
   */
  void remove_lambda_from_hessian(const double lambda);
  /*!
   * @brief Compute gain factor used in LM algorithm
   * @param cost_before_update least square cost before update vertices using delta_x_ in solver_ptr
   * @param lambda Lambda value
   * @return computed gain factor
   */
  double compute_gain_factor(const double cost_before_update, const double lambda);
private:
  /// base class pointer to solver
  std::shared_ptr<BaseSolver> solver_ptr_;
  /// Configuration for LM algorithm
  LevenbergMarguardtConfig config_;
};
}