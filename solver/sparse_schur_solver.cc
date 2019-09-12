#include "sparse_schur_solver.hh"
#include <iostream>
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
  int complement_params_size = 0;
  //First allocate start index for parameters in unmarginalized vertices
  for(auto it = problem_ptr_->vertex_id_map_vertex_ptr_.begin(); it!=problem_ptr_->vertex_id_map_vertex_ptr_.end(); it++){
    if(marginalized_vertices_.find(it->first) != marginalized_vertices_.end()){
      //a vertex need to be marginalized
      marginalized_params_size_+=it->second->minimal_dimension();
    } else{
      //a vertex not in the marginalized set
      vertex_id_map_start_index_[it->first] = complement_params_size;
      complement_params_size+=it->second->minimal_dimension();
    }
  }
  //Then allocate start index for parameters in vertices that need to be marginalized
  problem_params_size_ = complement_params_size;
  for(auto it = marginalized_vertices_.begin(); it!=marginalized_vertices_.end(); it++){
    if(problem_ptr_->vertex_id_map_vertex_ptr_.find(*it) != problem_ptr_->vertex_id_map_vertex_ptr_.end()){
      //Check if the vertex want to be marginalized is in the problem
      vertex_id_map_start_index_[*it] = problem_params_size_;
      problem_params_size_ += problem_ptr_->vertex_id_map_vertex_ptr_[*it]->minimal_dimension();
    }
  }
}

void SparseSchurSolver::compute_schur_complements() {
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
            int dim_compl = vertex_ptr_compl->minimal_dimension();
            shur_complement.second.emplace_back(start_index_compl, dim_compl);
          }
        }
      }
      schur_complements_.push_back(shur_complement);
    }
  }
}

void SparseSchurSolver::solve_delta_x()
{
  //New H_unmarg and g_unmarg
  int compl_size = problem_params_size_ - marginalized_params_size_;
  Eigen::MatrixXd hessian_compl = hessian_.block(0,0,compl_size, compl_size);
  Eigen::VectorXd g_compl = g_.segment(0,compl_size);
  Eigen::VectorXd g_marg = g_.segment(compl_size, marginalized_params_size_);

  //First perform marginalization
  for(auto& schur_complement : schur_complements_){
    Parameters marginalized = schur_complement.first;
    std::vector<Parameters> complements = schur_complement.second;
    Eigen::MatrixXd hessian_marg = hessian_.block(marginalized.start_index_, marginalized.start_index_,
                                                       marginalized.minimal_dim_, marginalized.minimal_dim_);
    Eigen::MatrixXd hessian_marg_inv = hessian_marg.inverse();
    //Perform marginalization to corresponding vertices
    for(int i = 0; i < complements.size(); i++){
      int start_index_i = complements[i].start_index_;
      int dim_i = complements[i].minimal_dim_;
      Eigen::MatrixXd block_i_mult_marg_inv = hessian_.block(start_index_i, marginalized.start_index_, dim_i, marginalized.minimal_dim_) * hessian_marg_inv;
      //Complement in gradient vector
      g_compl.segment(start_index_i, dim_i) -= block_i_mult_marg_inv * g_.segment(marginalized.start_index_, marginalized.minimal_dim_);
      //Complement in hessian matrix
      for(int j = i; j < complements.size(); j++){
        int start_index_j = complements[j].start_index_;
        int dim_j = complements[j].minimal_dim_;
        Eigen::MatrixXd complement_matrix = block_i_mult_marg_inv * hessian_.block(marginalized.start_index_, start_index_j, marginalized.minimal_dim_, dim_j);
        hessian_compl.block(start_index_i, start_index_j, dim_i, dim_j) -= complement_matrix;
        if(i!=j){
          hessian_compl.block(start_index_j, start_index_i, dim_j, dim_i) -= complement_matrix.transpose();
        }
      }
    }
  }
  //Solve complement unmarginalized parameters
  delta_x_.segment(0,compl_size) = hessian_compl.inverse() * g_compl;
  //Then solve for marginalized parameters
  g_marg -= hessian_.block(compl_size, 0, marginalized_params_size_, compl_size) * delta_x_.segment(0, compl_size);
  for(auto& schur_complement : schur_complements_) {
    Parameters marginalized = schur_complement.first;
    Eigen::MatrixXd hessian_marg_inv = hessian_.block(marginalized.start_index_, marginalized.start_index_,
                                                      marginalized.minimal_dim_, marginalized.minimal_dim_).inverse();
    delta_x_.segment(marginalized.start_index_, marginalized.minimal_dim_) = hessian_marg_inv * g_marg.segment(marginalized.start_index_-compl_size, marginalized.minimal_dim_);
  }
}



void SparseSchurSolver::set_marginalized_vertices(const std::set<IDType> &marginalized_vertices)
{
  marginalized_vertices_ = marginalized_vertices;
}

}