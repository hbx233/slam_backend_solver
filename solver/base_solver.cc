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

void BaseSolver::build_solve_structure() {
  // Make the Hessian Matrix and g vector in Optimization Normal Equation
  //Normal Equation: J.transpose() * J * delta_x = - J.transpose() * error
  // Initial Hessian Matrix and b vector are all zeros
  // Accumulate value on Hessian Matrix and b vector
  hessian_ = Eigen::MatrixXd::Zero(problem_params_size_, problem_params_size_);
  g_ = Eigen::VectorXd::Zero(problem_params_size_);

  //For Every Edge in the Problem
  for (auto &edge: problem_ptr_->edge_id_map_edge_ptr_) {
    //Compute errors and error function's jacobians
    edge.second->compute_errors();
    edge.second->compute_jacobians();
    //Row Block in Hessian
    for (size_t i = 0; i < edge.second->num_vertices(); i++) {
      std::shared_ptr<Vertex> v_ptr_i = edge.second->get_vertex_interface(i);
      if (v_ptr_i->is_fixed()) {
        continue;    // Hessian block is zero for fixed vertex
      }

      Eigen::MatrixXd& jacobian_i = edge.second->jacobians_[i];
      int start_index_i = vertex_id_map_start_index_[v_ptr_i->id()];
      int dim_i = v_ptr_i->minimal_dimension();
      Eigen::MatrixXd JtW = jacobian_i.transpose() * edge.second->information();
      //Column Block in Hessian
      //Start from i to exploit the symmetry in Hessian Matrix
      for (size_t j = i; j < edge.second->num_vertices(); j++) {
        std::shared_ptr<Vertex> v_ptr_j = edge.second->get_vertex_interface(j);

        if (v_ptr_j->is_fixed()){
          continue;
        }

        Eigen::MatrixXd& jacobian_j = edge.second->jacobians_[j];
        int start_index_j = vertex_id_map_start_index_[v_ptr_j->id()];
        int dim_j = v_ptr_j->minimal_dimension();
        Eigen::MatrixXd hessian = JtW * jacobian_j;
        hessian_.block(start_index_i,start_index_j,dim_i,dim_j).noalias() += hessian;
        //Symmetry block cross the diagonal line
        if (j != i) {
          hessian_.block(start_index_j,start_index_i,dim_j,dim_i).noalias() += hessian.transpose();
        }
      }
      g_.segment(start_index_i, dim_i).noalias() -= JtW * edge.second->errors_;
    }
  }
  delta_x_ = Eigen::VectorXd::Zero(problem_params_size_);  // initial delta_x = 0_n;
}


}

