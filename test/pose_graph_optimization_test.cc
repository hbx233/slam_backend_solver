//
// Created by hbx on 2019/9/14.
//
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
#include "solver/pose3d_vertex.hh"
#include "solver/relative_pose3d_edge.hh"
#include "solver/plain_solver.hh"
#include "solver/sparse_cholesky_solver.hh"
#include "solver/problem.hh"
#include "solver/minimizer_levenberg_marquardt.hh"

using namespace SLAMSolver;

//Generate single pose chain with loop closure
void get_simulation_data(std::vector<Sophus::SE3>& ground_truth, std::vector<Sophus::SE3>& relative_pose_obs){
  int pose_nums = 100;

  double radius = 8;

  ground_truth.clear();
  relative_pose_obs.clear();

  std::default_random_engine generator;
  std::normal_distribution<double> noise_rotation_pdf(0., 1. / 100);
  std::normal_distribution<double> noise_translation_pdf(0., 1/2);

  for (int n = 0; n <= pose_nums; ++n) {
    double theta = n * 2 * M_PI / pose_nums;

    Sophus::SO3 R(0,0,theta);
    Eigen::Vector3d t = Eigen::Vector3d(radius * cos(theta) - radius, radius * sin(theta), 1 * sin(2 * theta));
    ground_truth.push_back(Sophus::SE3{R,t});
    if(n >=1 ){
      Sophus::SE3 relative_pose = ground_truth[n-1].inverse() * ground_truth[n];
      Eigen::Vector3d noise_rotation(0,0,0);
      Eigen::Vector3d noise_translation(0,0,0);
      noise_rotation(0) += noise_rotation_pdf(generator);
      noise_rotation(1) += noise_rotation_pdf(generator);
      noise_rotation(2) += noise_rotation_pdf(generator);
      noise_translation(0) += noise_translation_pdf(generator);
      noise_translation(1) += noise_translation_pdf(generator);
      noise_translation(2) += noise_translation_pdf(generator);
      Sophus::SE3 noise_SE3(Sophus::SO3::exp(noise_rotation), noise_translation);
      relative_pose_obs.push_back(noise_SE3 * relative_pose);
    }
  }
}

int main(){
  std::vector<Sophus::SE3> poses_gt;
  std::vector<Sophus::SE3> relative_pose_obs;
  get_simulation_data(poses_gt, relative_pose_obs);

  //simply accumulate
  std::vector<Sophus::SE3> poses_accumulated;
  poses_accumulated.push_back(Sophus::SE3());
  for(int i = 0; i < relative_pose_obs.size(); i++){
    poses_accumulated.push_back(poses_accumulated[i] * relative_pose_obs[i]);
  }

  std::shared_ptr<Problem> problem_ptr = std::make_shared<Problem>();

  //Create 3D pose vertices, use accumulated pose as initial parameter
  std::vector<std::shared_ptr<Pose3DVertex>> pose_vertices;
  std::vector<unsigned long> solve_order;
  for(int i = 0; i < poses_accumulated.size() - 1; i++){
    std::shared_ptr<Pose3DVertex> pose_ptr = std::make_shared<Pose3DVertex>(true);
    pose_ptr->set_parameters(poses_accumulated[i]);
    if(i == 0){
      pose_ptr->set_fixed();
    }
    solve_order.push_back(pose_ptr->id());
    pose_vertices.push_back(pose_ptr);
    problem_ptr->add_vertex(pose_ptr);
  }
  std::cout<<pose_vertices.size()<<std::endl;
  Eigen::Matrix<double,6,6> info = Eigen::Matrix<double,6,6>::Zero();
  info(0,0) = 4;
  info(1,1) = 4;
  info(2,2) = 4;
  info(3,3) = 10000;
  info(4,4) = 10000;
  info(5,5) = 10000;
  //Create edge between two neighboring poses
  for(int i = 0; i < relative_pose_obs.size(); i++){
    std::shared_ptr<RelativePose3DEdge> edge_ptr = std::make_shared<RelativePose3DEdge>(relative_pose_obs[i]);
    edge_ptr->set_vertex<0>(pose_vertices[i]);
    if(i != relative_pose_obs.size() - 1){
      edge_ptr->set_vertex<1>(pose_vertices[i+1]);
    } else{
      //Loop closure constraints
      std::cout<<"Loop Closure"<<std::endl;
      edge_ptr->set_vertex<1>(pose_vertices[0]);
    }
    edge_ptr->set_information(info);
    problem_ptr->add_edge(edge_ptr);
  }

  std::cout<<"*"<<std::endl;
  //Create Problem Solver
  //std::shared_ptr<PlainSolver> solver_ptr = std::make_shared<PlainSolver>(problem_ptr);

  std::shared_ptr<SparseCholeskySolver> solver_chol_ptr = std::make_shared<SparseCholeskySolver>(problem_ptr);
  solver_chol_ptr->set_solve_order(solve_order);
#if 0
  solver_chol_ptr->set_solve_order(solve_order);
  solver_chol_ptr->compute_vertices_index();
  solver_chol_ptr->build_solve_structure();
  solver_chol_ptr->hessian_ += Eigen::MatrixXd::Identity(hessian_n, hessian_n);

  Eigen::LLT<Eigen::MatrixXd> llt_of_hessian(solver_chol_ptr->hessian_);
  Eigen::MatrixXd U = llt_of_hessian.matrixL();
  std::cout<<U<<std::endl;
  solver_chol_ptr->solve_delta_x();
  std::cout<<"*"<<std::endl;
  std::cout<<"delta_x from cholesky solver: "<<std::endl;
  std::cout<<solver_chol_ptr->delta_x_<<std::endl;
#endif
  //Create LM minimizer
#if 1
  LevenbergMarguardtConfig lm_config;
  MinimizerLevenbergMarquardt lm_minimizer(solver_chol_ptr, lm_config);
  lm_minimizer.minimize(1);

  for(int i = 0; i < poses_gt.size(); i++){
    std::cout<<"Number: "<<i<<"th pose in Pose Chain"<<std::endl;
    std::cout<<"=====Ground truth Pose====="<<std::endl;
    std::cout<<poses_gt[i].matrix()<<std::endl;
    std::cout<<"=====Accumulated Pose====="<<std::endl;
    std::cout<<poses_accumulated[i].matrix()<<std::endl;
    std::cout<<"=====Optimized Pose====="<<std::endl;
    std::cout<<pose_vertices[i]->parameters().matrix()<<std::endl;
  }
#endif
  return 0;
}