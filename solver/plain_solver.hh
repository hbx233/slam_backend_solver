//
// Created by hbx on 2019/9/9.
//
#include "solver/base_solver.hh"
#include "solver/problem.hh"
namespace SLAMSolver {
class PlainSolver : public BaseSolver {
public:
  PlainSolver(std::shared_ptr<Problem> problem_ptr);
  void compute_vertices_index() override;
  void build_solve_structure() override;
  void solve_delta_x() override;

private:
  int problem_variable_size_;
};
}