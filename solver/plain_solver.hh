//
// Created by hbx on 2019/9/9.
//
#include "solver/base_solver.hh"
#include "solver/problem.hh"
namespace SLAMSolver {
/*!
 * A plain Solver only use matrix inversion to solve normal equation
 */
class PlainSolver : public BaseSolver {
public:
  /// @brief Constructor of a plain solver
  PlainSolver(std::shared_ptr<Problem> problem_ptr);
  /// @brief Just use the ID order to compute the index
  void compute_vertices_index() override;
  /// @brief Just use matrix inversion to solve normal equation
  void solve_delta_x() override;
};
}