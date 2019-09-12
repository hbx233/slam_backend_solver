//
// Created by hbx on 2019/9/9.
//

#include "plain_solver.hh"
#include <iostream>
namespace SLAMSolver{
PlainSolver::PlainSolver(std::shared_ptr<Problem> problem_ptr)
: BaseSolver(problem_ptr)
{
}

void PlainSolver::compute_vertices_index(){
  //Default is using ID as relative order to construct whole optimization state vector
  //start_index start from zero
  problem_params_size_ = 0;
  for(auto it = problem_ptr_->vertex_id_map_vertex_ptr_.begin(); it!=problem_ptr_->vertex_id_map_vertex_ptr_.end(); it++){
    vertex_id_map_start_index_[it->first] = problem_params_size_;
    problem_params_size_ += it->second->minimal_dimension();
  }
}

void PlainSolver::solve_delta_x() {
  delta_x_ = hessian_.inverse() * g_;
}

}