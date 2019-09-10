#pragma once

#include "solver/base_solver.hh"
#include <set>
namespace SLAMSolver{

struct Parameters{
  Parameters(int start_index, int minimal_dim)
  : start_index_(start_index), minimal_dim_(minimal_dim)
  {}
  int start_index_;
  int minimal_dim_;
};

//Assume there's no edge between any pair of vertices that want to be marginalized
class SparseSchurSolver : public BaseSolver{
public:
  SparseSchurSolver(std::shared_ptr<Problem> problem_ptr);
  void compute_vertices_index() override;
  void build_solve_structure() override;
  void solve_delta_x() override;
  void set_marginalized_vertices(const std::set<IDType>& marginalized_vertices);
private:
  void compute_marginalization_operations();
  void perform_block_schur_complement();
  int problem_variable_size_;
  int marginalized_variable_size_;
  std::set<IDType> marginalized_vertices_;
  std::vector<std::pair<Parameters, std::vector<Parameters>>> shur_complement_blocks_;
};
}
