#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "solver/problem.hh"
#include "utils/tic_toc.h"

namespace SLAMSolver {

Problem::Problem()
{
}

Problem::~Problem() {}

bool Problem::add_vertex(std::shared_ptr<Vertex> vertex) {
  if (vertex_id_map_vertex_ptr_.find(vertex->id()) != vertex_id_map_vertex_ptr_.end()) {
    return false;
  } else {
    vertex_id_map_vertex_ptr_[vertex->id()] = vertex;
    return true;
  }
}


bool Problem::add_edge(std::shared_ptr<Edge> edge) {
  //Add Edge Interface pointer to hash map with its ID as key 
  if (edge_id_map_edge_ptr_.find(edge->id()) != edge_id_map_edge_ptr_.end()) {
    return false;
  } else {
    edge_id_map_edge_ptr_[edge->id()] = edge;
    for (int i = 0; i < edge->num_vertices(); i++){
      auto vertex_interface_ptr = edge->get_vertex_interface(i);
      add_vertex(vertex_interface_ptr);
    }
    return true;
  }
}

void Problem::compute_start_index(){
  //Default is using ID as relative order to construct whole optimization state vector
  //start_index start from zero
  for(auto it = vertex_id_map_vertex_ptr_.begin(); it!=vertex_id_map_vertex_ptr_.end(); it++){
    vertex_id_map_start_index_[it->first] = problem_variable_size_;
    problem_variable_size_ += it->second->minimal_dimension();
  }
}

void Problem::build_gauss_newton() {
  // Make the Hessian Matrix and b vector in Optimization Normal Equation
  //Normal Equation: J.transpose() * J * delta_x = - J.transpose() * error
  // Initial Hessian Matrix and b vector are all zeros
  // Accumulate value on Hessian Matrix and b vector 
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(problem_variable_size_, problem_variable_size_));
  Eigen::VectorXd b(Eigen::VectorXd::Zero(problem_variable_size_));
  
  //For Every Edge in the Problem
  for (auto &edge: edge_id_map_edge_ptr_) {
    //Compute errors and error function's jacobians
    edge.second->compute_errors();
    edge.second->compute_jacobians();
    //Row Block in Hessian
    //std::cout<<"Edge: "<<edge.second->id()<<"=========="<<std::endl;
    for (size_t i = 0; i < edge.second->num_vertices(); i++) {
      std::shared_ptr<Vertex> v_ptr_i = edge.second->get_vertex_interface(i);
      if (v_ptr_i->is_fixed()) {
        continue;    // Hessian block is zero for fixed vertex
      }

      Eigen::MatrixXd& jacobian_i = edge.second->jacobians_[i];
      //std::cout<<"Jacobian of: "<<i<<"th"<<" vertex"<<std::endl;
      //std::cout<<jacobian_i<<std::endl;
      int start_index_i = vertex_id_map_start_index_[v_ptr_i->id()];
      int dim_i = v_ptr_i->minimal_dimension();

      //std::cout<<edge.second->information()<<std::endl;
      Eigen::MatrixXd JtW = jacobian_i.transpose() * edge.second->information();
	
      //Column Block in Hessian
      //Start from i to exploit the symmetry in Hessian Matrix
      for (size_t j = i; j < edge.second->num_vertices(); j++) {
        std::shared_ptr<Vertex> v_ptr_j = edge.second->get_vertex_interface(j);

        if (v_ptr_j->is_fixed()){
          continue;
        }

        Eigen::MatrixXd& jacobian_j = edge.second->jacobians_[j];
        //std::cout<<"Jacobian of: "<<j<<"th"<<" vertex"<<std::endl;
        //std::cout<<jacobian_j<<std::endl;
        int start_index_j = vertex_id_map_start_index_[v_ptr_j->id()];
        int dim_j = v_ptr_j->minimal_dimension();

        Eigen::MatrixXd hessian = JtW * jacobian_j;
        //std::cout<<"Hessian to be add: "<<hessian<<std::endl;
        H.block(start_index_i,start_index_j,dim_i,dim_j).noalias() += hessian;
	       //Symmetry block cross the diagonal line
        if (j != i) {
          H.block(start_index_j,start_index_i,dim_j,dim_i).noalias() += hessian.transpose();
        }
      }
      //std::cout<<"Edge error: "<<edge.second->errors_<<std::endl;
      //std::cout<<"b vector to be add: "<<JtW * edge.second->errors_<<std::endl;
      b.segment(start_index_i, dim_i).noalias() += JtW * edge.second->errors_;
    }
  }
  Hessian_ = H;
  b_ = b;


//    Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
//    std::cout << svd.singularValues() <<std::endl;

//  if (err_prior_.rows() > 0) {
//      b_prior_ -= H_prior_ * delta_x_.head(ordering_poses_);   // update the error_prior
//  }
//  Hessian_.topLeftCorner(ordering_poses_, ordering_poses_) += H_prior_;
//  b_.head(ordering_poses_) += b_prior_;

  delta_x_ = Eigen::VectorXd::Zero(problem_variable_size_);  // initial delta_x = 0_n;
}


void Problem::solve_normal_equation() {
  delta_x_ = - Hessian_.inverse() * b_;
}

double Problem::compute_initial_lambda(const Eigen::MatrixXd &Hessian, const double tau) {
  //Check Hessian Matrix
  if(Hessian.cols() != Hessian.rows()){
    std::cout<<"[ERROR] Hessian Matrix is not Square"<<std::endl;
    exit(0);
  }
  double max_hessian_diag = 0;
  for(int i = 0; i < Hessian.cols(); i++){
    max_hessian_diag = std::max(max_hessian_diag, Hessian(i,i));
  }
  return tau * max_hessian_diag;
}

void Problem::add_lambda_to_hessian(const double lambda) {
  //Check Hessian Matrix
  if(Hessian_.cols() != Hessian_.rows()){
    std::cout<<"[ERROR] Hessian Matrix is not Square"<<std::endl;
    exit(0);
  }
  for(int i = 0; i < Hessian_.cols(); i++) {
    Hessian_(i,i) += lambda;
  }
}

void Problem::remove_lambda_from_hessian(const double lambda) {
  //Check Hessian Matrix
  if(Hessian_.cols() != Hessian_.rows()){
    std::cout<<"[ERROR] Hessian Matrix is not Square"<<std::endl;
    exit(0);
  }
  for(int i = 0; i < Hessian_.cols(); i++) {
    Hessian_(i,i) -= lambda;
  }
}

void Problem::update_vertices() {
  for(auto it = vertex_id_map_vertex_ptr_.begin(); it!=vertex_id_map_vertex_ptr_.end(); it++){
    auto v_ptr = it->second;
    int start_index = vertex_id_map_start_index_[it->first];
    int dim = v_ptr->minimal_dimension();
    //Fetch this vertex's corresponding update vector from delta_x_
    Eigen::VectorXd delta = delta_x_.segment(start_index,dim);
    //Update vertex through plus interface
    v_ptr->plus(delta);
  }
}

void Problem::rollback_vertices() {
  for(auto it = vertex_id_map_vertex_ptr_.begin(); it!=vertex_id_map_vertex_ptr_.end(); it++){
    auto v_ptr = it->second;
    int start_index = vertex_id_map_start_index_[it->first];
    int dim = v_ptr->minimal_dimension();
    Eigen::VectorXd delta = delta_x_.segment(start_index,dim);
    v_ptr->plus(-delta);
  }
}

double Problem::compute_cost() {
  //Compute total cost
  double total_cost = 0;
  for (auto &edge: edge_id_map_edge_ptr_) {
    //Compute errors and error function's jacobians
    edge.second->compute_errors();
    total_cost += edge.second->chi2();
  }
  return total_cost;
}

double Problem::compute_gain_factor(const double cost_before_update, const double lambda) {
  //compute cost after update
  double cost_after_update = compute_cost();
  //compute Gauss-Newton approximated cost
  double gn_cost = delta_x_.transpose() * (lambda * delta_x_ - b_) + 0.001;
  return (cost_before_update - cost_after_update) / gn_cost;
}

bool Problem::solve(int max_iterations) {
  int iter = 0; //optimization iteration
  double v = 2;
  //Gauss-Newton, build Hessian Matrix and b vector in Normal Equation
  compute_start_index();
  build_gauss_newton();
  //Compute Initial value of lambda in LM
  double lambda = compute_initial_lambda(Hessian_,0.000001);
  std::cout<<"Initial Lambda: "<<lambda<<std::endl;
  //Threshold for stopping criteria
  double th_1 = 0.00000001;
  double th_2 = 0.000001;
  bool found = b_.lpNorm<Eigen::Infinity>() < th_1;
  while(!found && iter < max_iterations){
    iter++;
    bool update_successful = false;
    int fail_count = 0;
    while(!update_successful){
      //1. Add lambda to Hessian and construct the normal equation
      std::cout<<"Current Lambda: "<<lambda<<std::endl;
      add_lambda_to_hessian(lambda);
      //2. Solve Linear system to compute delta_x_:
      // (Hessian + lambda * I) * delta_x = -b
      solve_normal_equation();
      std::cout<<"delta_x: "<<std::endl<<delta_x_<<std::endl;
      if(delta_x_.norm() <= th_2 || fail_count > 10){
        //Stop Criteria: delta_x_ is too small
        found = true;
        break;
      } else {
        //1. Cache total cost before update
        double cost_before_update = compute_cost();
        //2. Update all vertices using delta_x_
        update_vertices();
        //3. Compute gain factor
        double gain = compute_gain_factor(cost_before_update, lambda);
        if (gain > 0) {
          std::cout << "A Good Step" << std::endl;
          //A good update, cost goes down
          update_successful = true;
          fail_count = 0;
          //1. Check if satisfied the first stopping criteria
          if (b_.lpNorm<Eigen::Infinity>() < th_1) {
            found = true;
            continue;
          }
          //2. Compute new Hessian and b vector with Gauss-Newton
          build_gauss_newton();
          //3. Update Lambda
          double factor = 1 - std::pow(2 * gain - 1, 3);
          lambda = lambda * std::max(0.333, factor);
          v = 2;
        } else {
          std::cout << "Not a good step, roll back" << std::endl;
          //Not a good step
          update_successful = false;
          fail_count++;
          // 1.roll back to original parameter
          rollback_vertices();
          // 2.Remove previous lambda from Hessian Matrix
          remove_lambda_from_hessian(lambda);
          // 3. Update lambda and v;
          lambda *= v;
          v *= 2;
        }
      }
    }
  }
  return true;
}

}






