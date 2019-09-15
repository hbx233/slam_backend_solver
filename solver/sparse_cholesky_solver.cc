//
// Created by hbx on 2019/9/14.
//

#include "solver/sparse_cholesky_solver.hh"

namespace SLAMSolver{
SparseCholeskySolver::SparseCholeskySolver(std::shared_ptr<Problem> problem_ptr)
: BaseSolver(problem_ptr){}

void SparseCholeskySolver::compute_vertices_index() {
  problem_params_size_ = 0;
  for(auto& order : solve_order_){
    block_dim_ = problem_ptr_->vertex_id_map_vertex_ptr_[order]->minimal_dimension();
    vertex_id_map_start_index_[order] = problem_params_size_;
    problem_params_size_ += problem_ptr_->vertex_id_map_vertex_ptr_[order]->minimal_dimension();
  }
}

void SparseCholeskySolver::set_solve_order(const std::vector<IDType>& solve_order) {
  solve_order_ = solve_order;
}

void SparseCholeskySolver::solve_delta_x() {
  //number of blocks
  int n = solve_order_.size();
  std::vector<std::vector<Eigen::MatrixXd>> cholesky_L(n, std::vector<Eigen::MatrixXd>(n));
  //sparse cholesky block decomposition
  for(int i = 1; i <n; i++){
    Eigen::MatrixXd h_block = hessian_.block(i*block_dim_, i*block_dim_,block_dim_,block_dim_);
    if(i == 1){
      Eigen::LLT<Eigen::MatrixXd> llt_of_h_block(h_block);
      cholesky_L[i][i] = llt_of_h_block.matrixL();
    } else{
      cholesky_L[i][i-1] = hessian_.block(i*block_dim_, (i-1)*block_dim_, block_dim_, block_dim_) * cholesky_L[i-1][i-1].transpose().inverse();
      Eigen::MatrixXd new_block = h_block - cholesky_L[i][i-1] * cholesky_L[i][i-1].transpose();
      Eigen::LLT<Eigen::MatrixXd> llt_of_h_block(new_block);
      cholesky_L[i][i] = llt_of_h_block.matrixL();
    }
  }
  //back-ward pass to solve U * c = g
  Eigen::VectorXd c;
  c.resize(problem_params_size_);
  for(int i = 1; i <n; i++){
    if(i == 1){
      c.segment(i*block_dim_,block_dim_) = cholesky_L[i][i].inverse() * g_.segment(i*block_dim_, block_dim_);
    } else{
      c.segment(i*block_dim_, block_dim_) = cholesky_L[i][i].inverse() * (g_.segment(i*block_dim_, block_dim_) - cholesky_L[i][i-1]*c.segment((i-1)*block_dim_,block_dim_));
    }
  }
  //forward-pass, solve delta_x_
  for(int i = n-1; i > 0; i--){
    if(i == n-1){
      delta_x_.segment(i*block_dim_,block_dim_) = cholesky_L[i][i].transpose().inverse() * c.segment(i*block_dim_,block_dim_);
    } else{
      delta_x_.segment(i*block_dim_,block_dim_) = cholesky_L[i][i].transpose().inverse() * (c.segment(i*block_dim_,block_dim_) - cholesky_L[i+1][i].transpose() * delta_x_.segment((i+1)*block_dim_,block_dim_));
    }
  }
}


}

