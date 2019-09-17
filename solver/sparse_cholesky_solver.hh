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
  /// @brief Compute Vertices index with provided solve order
  /// @note Different solve order will make problem with different solving efficiency
  void compute_vertices_index() override;
  /// @brief Solve delta_x_ using sparse cholesky decomposition
  void solve_delta_x() override;
  /// @brief Build the solve structure and the sparse structure of hessian matrix
  void build_solve_structure() override;
  /*!
   * @brief Set the order of vertices to solve more efficiently
   * @param solve_order The vertices order in solving
   */
  void set_solve_order(const std::vector<IDType>& solve_order);
private:
  int block_dim_;
  std::vector<IDType> solve_order_;
  std::map<IDType, int> vertex_id_map_order_;
  std::vector<std::vector<Block>> hessian_block_;
};
}
