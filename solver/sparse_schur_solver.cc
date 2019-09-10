//
// Created by hbx on 2019/9/9.
//

#include "sparse_schur_solver.hh"

namespace SLAMSolver{
SparseSchurSolver::SparseSchurSolver(std::shared_ptr<Problem> problem_ptr)
: BaseSolver(problem_ptr)
{
}

void SparseSchurSolver::compute_vertices_index()
{
  //Order the problem as
  // complement block (unmarginalized) in top-left big diagonal block
  // marginalized block in bottom-right big diagonal block
  int unmarginalized_variable_size = 0;
  //First allocate start index for parameters in unmarginalized vertices
  for(auto it = problem_ptr_->vertex_id_map_vertex_ptr_.begin(); it!=problem_ptr_->vertex_id_map_vertex_ptr_.end(); it++){
    if(marginalized_vertices_.find(it->first) != marginalized_vertices_.end()){
      //a vertex need to be marginalized
      marginalized_variable_size_+=it->second->minimal_dimension();
    } else{
      //a vertex not in the marginalized set
      vertex_id_map_start_index_[it->first] = unmarginalized_variable_size;
      unmarginalized_variable_size+=it->second->minimal_dimension();
    }
  }
  //Then allocate start index for parameters in vertices that need to be marginalized
  problem_variable_size_ = marginalized_variable_size_;
  for(auto it = marginalized_vertices_.begin(); it!=marginalized_vertices_.end(); it++){
    if(problem_ptr_->vertex_id_map_vertex_ptr_.find(*it) != problem_ptr_->vertex_id_map_vertex_ptr_.end()){
      //Check if the vertex want to be marginalized is in the problem
      vertex_id_map_start_index_[*it] = problem_variable_size_;
      problem_variable_size_ += problem_ptr_->vertex_id_map_vertex_ptr_[*it]->minimal_dimension();
    }
  }
}

void SparseSchurSolver::build_solve_structure()
{
  // Make the Hessian Matrix and g vector in Optimization Normal Equation
  //Normal Equation: J.transpose() * J * delta_x = - J.transpose() * error
  // Initial Hessian Matrix and b vector are all zeros
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

void SparseSchurSolver::compute_marginalization_operations() {
  //Loop through every Vertex want to be marginalized
  for(auto it = marginalized_vertices_.begin(); it != marginalized_vertices_.end(); it++){
    IDType vertex_id_marg = *it;
    if(problem_ptr_->vertex_id_map_edges_ptr_.find(vertex_id_marg) != problem_ptr_->vertex_id_map_edges_ptr_.end()){
      //Exist in the problem
      auto vertex_ptr_marg = problem_ptr_->vertex_id_map_vertex_ptr_[vertex_id_marg];
      int dim_marg = vertex_ptr_marg->minimal_dimension();
      int start_index_marg = vertex_id_map_start_index_[vertex_id_marg];
      auto connected_edges = problem_ptr_->vertex_id_map_edges_ptr_.equal_range(*it);
      //Create shur complement block for the marginalized vertex
      Parameters marg_block(start_index_marg, dim_marg);
      std::pair<Parameters, std::vector<Parameters>> shur_complement(marg_block, std::vector<Parameters>());
      //For every connected edges
      for(auto it_edge = connected_edges.first; it_edge != connected_edges.second; it_edge++){
        //For every vertices that connected with the marginalized vertex through this edge
        auto edge_ptr = it_edge->second;
        for(int i = 0; i < edge_ptr->num_vertices(); i++){
          auto vertex_ptr_compl = edge_ptr->get_vertex_interface(i);
          //If vertex is not the marginalized vertex, add to marginalization complement
          if(vertex_ptr_compl->id() != vertex_id_marg){
            int start_index_compl = vertex_id_map_start_index_[vertex_ptr_compl->id()];
            int dim_compl = vertex_ptr_compl->id();
            shur_complement.second.emplace_back(start_index_compl, dim_compl);
          }
        }
      }
    }
  }
}

void SparseSchurSolver::solve_delta_x()
{
  //New H_unmarg and g_unmarg
  int unmarg_size = problem_variable_size_ - marginalized_variable_size_;
  Eigen::MatrixXd hessian_unmarg = hessian_.block(0,0,unmarg_size, unmarg_size);
  Eigen::VectorXd g_unmarg = g_.segment(0,unmarg_size);
  Eigen::VectorXd g_marg = g_.segment(unmarg_size, marginalized_variable_size_);
  //First perform marginalization
  for(auto& shur_complement_block : shur_complement_blocks_){
    Parameters marginalized = shur_complement_block.first;
    std::vector<Parameters> complements = shur_complement_block.second;
    Eigen::MatrixXd hessian_marg_inv = hessian_.block(marginalized.start_index_, marginalized.start_index_,
                                                       marginalized.minimal_dim_, marginalized.minimal_dim_).inverse();
    for(int i = 0; i < complements.size(); i++){
      Eigen::MatrixXd
      for(int j = i; j < complements.size(); j++){
        Eigen::MatrixXd complement_matrix = ;
      }
    }
  }
  //Solve unmarginalized parameters
  Eigen::VectorXd delta_x_unmarg = hessian_unmarg.inverse() * g_unmarg;
  //Then solve for marginalized parameters
  
}



void SparseSchurSolver::set_marginalized_vertices(const std::set<IDType> &marginalized_vertices)
{
  marginalized_vertices_ = marginalized_vertices;
}

}