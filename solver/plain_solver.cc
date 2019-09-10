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
  for(auto it = problem_ptr_->vertex_id_map_vertex_ptr_.begin(); it!=problem_ptr_->vertex_id_map_vertex_ptr_.end(); it++){
    vertex_id_map_start_index_[it->first] = problem_variable_size_;
    problem_variable_size_ += it->second->minimal_dimension();
  }
}

void PlainSolver::build_solve_structure() {
  // Make the Hessian Matrix and g vector in Optimization Normal Equation
  //Normal Equation: J.transpose() * J * delta_x = - J.transpose() * error
  // H * delta_x = g
  // Initial Hessian Matrix and g vector are all zeros
  // Accumulate value on Hessian Matrix and b vector
  hessian_ = Eigen::MatrixXd::Zero(problem_variable_size_, problem_variable_size_);
  g_ = Eigen::VectorXd::Zero(problem_variable_size_);

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
#if 0
        std::cout<<"Hessian Contribution from: Edge "<<edge.first<<std::endl;
        std::cout<<hessian<<std::endl;
        std::cout<<"Jacobian: "<<std::endl;
        std::cout<<jacobian_i<<std::endl;
        std::cout<<"Information: "<<std::endl;
        std::cout<<edge.second->information()<<std::endl;
#endif
        hessian_.block(start_index_i,start_index_j,dim_i,dim_j).noalias() += hessian;
        //Symmetry block cross the diagonal line
        if (j != i) {
          hessian_.block(start_index_j,start_index_i,dim_j,dim_i).noalias() += hessian.transpose();
        }
      }

      g_.segment(start_index_i, dim_i).noalias() -= JtW * edge.second->errors_;
    }
  }
  delta_x_ = Eigen::VectorXd::Zero(problem_variable_size_);  // initial delta_x = 0_n;
}

void PlainSolver::solve_delta_x() {
  delta_x_ = hessian_.inverse() * g_;
}
}