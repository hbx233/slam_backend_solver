#pragma once
#include "solver/eigen_types.h"
#include "solver/problem.hh"
namespace SLAMSolver{
class BaseSolver{
public:
  using IDType = Problem::IDType;
  BaseSolver(std::shared_ptr<Problem> problem_ptr);
  void update_vertices_with_delta_x();
  void rollback_verteices();
  void modify_hessian(std::function<void(const Eigen::MatrixXd&)> modify_func);
  double compute_total_cost();

  virtual void compute_vertices_index() = 0;
  virtual void build_solve_structure() = 0;
  virtual void solve_delta_x() = 0;

public:
  Eigen::MatrixXd hessian_;
  Eigen::VectorXd delta_x_;
  Eigen::VectorXd g_;
protected:
  // Hash map from Vertex's ID to Vertex's start index in Problem's
  std::map<IDType, int> vertex_id_map_start_index_;
  std::shared_ptr<Problem> problem_ptr_;
};
}