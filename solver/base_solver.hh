#pragma once
#include "solver/eigen_types.h"
#include "solver/problem.hh"
namespace SLAMSolver{
/*!
 * @brief Base Class for all kinds of solvers
 */
class BaseSolver{
public:
  using IDType = Problem::IDType;
  /*!
   * @brief Constructor for BaseSolver
   * @param problem_ptr The problem need to be solved
   */
  BaseSolver(std::shared_ptr<Problem> problem_ptr);

  /*!
   * @brief Update all vertices' parameters using the current delta_x_
   */
  void update_vertices_with_delta_x();

  /*!
   * @brief Rollback all vertices' parameters for delta_x_ change
   */
  void rollback_verteices();

  /*!
   * @breif Compute current problem's total cost
   * @return total cost
   */
  double compute_total_cost();

  /*!
   * @brief Pure Virtual function for computing each vertex's parameter vector's start index
   *     Different solver need to implement their own strategy for ordering the vertics
   */
  virtual void compute_vertices_index() = 0;
  /*!
   * @brief Virtual function for building the normal equation given the vertices' position
   *    Compute hessian_ and g_
   *    Base class implemented a common version of building normal equation
   */
  virtual void build_solve_structure();
  /*!
   * @brief Pure Virtual function for solving the normal equation
   *        hessian_ * delta_x_ = g_
   */
  virtual void solve_delta_x() = 0;

public:
  Eigen::MatrixXd hessian_;
  Eigen::VectorXd delta_x_;
  Eigen::VectorXd g_;
protected:
  int problem_params_size_ = 0;
  // Hash map from Vertex's ID to Vertex's start index in Problem's
  std::map<IDType, int> vertex_id_map_start_index_;
  std::shared_ptr<Problem> problem_ptr_;
};
}