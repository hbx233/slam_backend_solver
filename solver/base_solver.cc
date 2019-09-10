#include "solver/base_solver.hh"

namespace SLAMSolver{
BaseSolver::BaseSolver(std::shared_ptr<Problem> problem_ptr)
:problem_ptr_(problem_ptr){}

void BaseSolver::update_vertices_with_delta_x() {
  for(auto it = problem_ptr_->vertex_id_map_vertex_ptr_.begin(); it!=problem_ptr_->vertex_id_map_vertex_ptr_.end(); it++){
    auto v_ptr = it->second;
    int start_index = vertex_id_map_start_index_[it->first];
    int dim = v_ptr->minimal_dimension();
    //Fetch this vertex's corresponding update vector from delta_x_
    Eigen::VectorXd delta = delta_x_.segment(start_index,dim);
    //Update vertex through plus interface
    v_ptr->plus(delta);
  }
}

void BaseSolver::rollback_verteices() {
  for(auto it = problem_ptr_->vertex_id_map_vertex_ptr_.begin(); it!=problem_ptr_->vertex_id_map_vertex_ptr_.end(); it++){
    auto v_ptr = it->second;
    int start_index = vertex_id_map_start_index_[it->first];
    int dim = v_ptr->minimal_dimension();
    Eigen::VectorXd delta = delta_x_.segment(start_index,dim);
    v_ptr->plus(-delta);
  }
}

double BaseSolver::compute_total_cost() {
  return problem_ptr_->compute_cost();
}



}

