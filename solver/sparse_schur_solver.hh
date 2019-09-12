#pragma once

#include "solver/base_solver.hh"
#include <set>
namespace SLAMSolver{

/// Parameter blocks in hessian matrix
struct Parameters{
  Parameters(int start_index, int minimal_dim)
  : start_index_(start_index), minimal_dim_(minimal_dim)
  {}
  int start_index_;
  int minimal_dim_;
};

/*!
 * @brief Solver using sparse schur decomposition to solve normal equation
 * Assume there's no edge between any pair of vertices that want to be marginalized
 */
class SparseSchurSolver : public BaseSolver{
public:
  /// @brief Constructor
  SparseSchurSolver(std::shared_ptr<Problem> problem_ptr);

  /// @brief Reorder the vertices, make all parameters of marginalized vertices to bottom-right
  void compute_vertices_index() override;
  /// @brief Solve normal equation using sparse schur complement
  void solve_delta_x() override;
  /*!
   * @brief Set vertices that want to marginalize
   * @param marginalized_vertices A set containing id of all vertices that want to marginalize
   */
  void set_marginalized_vertices(const std::set<IDType>& marginalized_vertices);
public:
  /// @brief Compute sparse schur complements that need to be performed on hessian matrix and g vector
  void compute_schur_complements();

  int marginalized_params_size_; // size of parameters that want to marginalize
  std::set<IDType> marginalized_vertices_; //a set contain id of vertices that want to marginalize
  /// represent the schur complement operation
  /// First in pair represent the parameters that need to be marginalized out
  /// Second in pair represent a set of parameters that need to complement when marginalizing out the pair.first
  std::vector<std::pair<Parameters, std::vector<Parameters>>> schur_complements_;
};
}
