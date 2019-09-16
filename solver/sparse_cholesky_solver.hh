#pragma once

#include "solver/base_solver.hh"
namespace SLAMSolver{

struct Block{
  Block(int r, int c, int dim_r, int dim_c, bool occupied)
  : r(r), c(c), dim_r(dim_r), dim_c(dim_c), occupied(occupied)
  {}
  Block() = default;
  ~Block() = default;
  int r = -1;
  int c = -1;
  int dim_r = -1;
  int dim_c = -1;
  bool occupied = false;
};

class SparseCholeskySolver : public BaseSolver{
public:
  /// @brief Constructor passes problem pointer to BaseSolver
  SparseCholeskySolver(std::shared_ptr<Problem> problem_ptr);
  /// @brief Default destructor
  ~SparseCholeskySolver() = default;
  void compute_vertices_index() override;
  void solve_delta_x() override;
  void build_solve_structure() override;
  void set_solve_order(const std::vector<IDType>& solve_order);
private:
  int block_dim_;
  std::vector<IDType> solve_order_;
  std::map<IDType, int> vertex_id_map_order_;
  std::vector<std::vector<Block>> hessian_block_;
};
}
