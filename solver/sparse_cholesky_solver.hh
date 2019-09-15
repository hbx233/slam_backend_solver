#pragma once

#include "solver/base_solver.hh"
namespace SLAMSolver{
class SparseCholeskySolver : public BaseSolver{
public:
  /// @brief Constructor passes problem pointer to BaseSolver
  SparseCholeskySolver(std::shared_ptr<Problem> problem_ptr);
  /// @brief Default destructor
  ~SparseCholeskySolver() = default;
  void compute_vertices_index() override;
  void solve_delta_x() override;
  void set_solve_order(const std::vector<IDType>& solve_order);
private:
  int block_dim_;
  std::vector<IDType> solve_order_;
};
}
