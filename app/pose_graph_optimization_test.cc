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

struct RelativePoseObservation{
  RelativePoseObservation(int pose0_idx, int pose1_idx, Sophus::SE3 pose0_SE3_pose1)
  : pose0_idx(pose0_idx), pose1_idx(pose1_idx), pose0_SE3_pose1(pose0_SE3_pose1){}
  int pose0_idx;
  int pose1_idx;
  Sophus::SE3 pose0_SE3_pose1;
};

Sophus::SE3 generate_SE3_noise(){
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::normal_distribution<double> noise_rotation_pdf(0., 1. / 100);
  std::normal_distribution<double> noise_x_pdf(0., 1/2);
  std::normal_distribution<double> noise_y_pdf(0., 1/2);
  std::normal_distribution<double> noise_z_pdf(0., 1/5);
  Eigen::Vector3d noise_rotation(0,0,0);
  Eigen::Vector3d noise_translation(0,0,0);
  noise_rotation(0) += noise_rotation_pdf(generator);
  noise_rotation(1) += noise_rotation_pdf(generator);
  noise_rotation(2) += noise_rotation_pdf(generator);
  noise_translation(0) += noise_x_pdf(generator);
  noise_translation(1) += noise_y_pdf(generator);
  noise_translation(2) += noise_z_pdf(generator);
  return Sophus::SE3(Sophus::SO3::exp(noise_rotation), noise_translation);
}

//Generate single pose chain with loop closure
void get_simulation_data(std::vector<Sophus::SE3>& ground_truth, std::vector<RelativePoseObservation>& relative_pose_obs){
  int pose_nums = 100;

  double radius = 15;

  ground_truth.clear();
  relative_pose_obs.clear();


  for (int n = 0; n <= pose_nums; ++n) {
    double theta = n * 2 * M_PI / pose_nums;

    Sophus::SO3 R(0,0,theta);
    Eigen::Vector3d t = Eigen::Vector3d(radius * cos(theta) - radius, radius * sin(theta), 1 * sin(2 * theta));
    ground_truth.push_back(Sophus::SE3{R,t});
  }
  for (int n = 0; n < ground_truth.size(); n++){
    if(n >=1 ){
      Sophus::SE3 relative_pose = ground_truth[n-1].inverse() * ground_truth[n];
      Sophus::SE3 noise_SE3 = generate_SE3_noise();
      RelativePoseObservation obs(n-1, n, noise_SE3 * relative_pose);
      relative_pose_obs.push_back(obs);
    }
  }
  for (int n = 0; n < ground_truth.size(); n++){
    if(n >= 3){
      Sophus::SE3 relative_pose = ground_truth[n-3].inverse() * ground_truth[n];
      Sophus::SE3 noise_SE3 = generate_SE3_noise();
      RelativePoseObservation obs(n-3, n, noise_SE3 * relative_pose);
      relative_pose_obs.push_back(obs);
    }
  }
  for (int n = 0; n < ground_truth.size(); n++){
    if(n >= 5){
      Sophus::SE3 relative_pose = ground_truth[n-5].inverse() * ground_truth[n];
      Sophus::SE3 noise_SE3 = generate_SE3_noise();
      RelativePoseObservation obs(n-5, n, noise_SE3 * relative_pose);
      relative_pose_obs.push_back(obs);
    }
  }
  for (int n = 0; n < ground_truth.size(); n++){
    if(n >= 9){
      Sophus::SE3 relative_pose = ground_truth[n-9].inverse() * ground_truth[n];
      Sophus::SE3 noise_SE3 = generate_SE3_noise();
      RelativePoseObservation obs(n-9, n, noise_SE3 * relative_pose);
      relative_pose_obs.push_back(obs);
    }
  }
}

int main(){
  std::vector<Sophus::SE3> poses_gt;
  std::vector<RelativePoseObservation> relative_pose_obs;
  get_simulation_data(poses_gt, relative_pose_obs);

  //simply accumulate
  std::vector<Sophus::SE3> poses_accumulated;
  poses_accumulated.push_back(Sophus::SE3());
  for(int i = 0; i < poses_gt.size()-1; i++){
    poses_accumulated.push_back(poses_accumulated[i] * relative_pose_obs[i].pose0_SE3_pose1);
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
  info(2,2) = 25;
  info(3,3) = 10000;
  info(4,4) = 10000;
  info(5,5) = 10000;
  //Create edge between two neighboring poses
  for(int i = 0; i < relative_pose_obs.size(); i++){
    std::shared_ptr<RelativePose3DEdge> edge_ptr = std::make_shared<RelativePose3DEdge>(relative_pose_obs[i].pose0_SE3_pose1);
    int pose0_idx = relative_pose_obs[i].pose0_idx;
    int pose1_idx = relative_pose_obs[i].pose1_idx;
    pose0_idx = pose0_idx == pose_vertices.size() ? 0 : pose0_idx;
    pose1_idx = pose1_idx == pose_vertices.size() ? 0 : pose1_idx;
    edge_ptr->set_vertex<0>(pose_vertices[pose0_idx]);
    edge_ptr->set_vertex<1>(pose_vertices[pose1_idx]);
    edge_ptr->set_information(info);
    problem_ptr->add_edge(edge_ptr);
  }

  //Create Problem Solver
  std::shared_ptr<PlainSolver> solver_ptr = std::make_shared<PlainSolver>(problem_ptr);

  std::shared_ptr<SparseCholeskySolver> solver_chol_ptr = std::make_shared<SparseCholeskySolver>(problem_ptr);
  solver_chol_ptr->set_solve_order(solve_order);
  //Create LM minimizer
  LevenbergMarguardtConfig lm_config;
  MinimizerLevenbergMarquardt lm_minimizer(solver_chol_ptr, lm_config);
  lm_minimizer.minimize(30);

  for(int i = 0; i < poses_gt.size()-1; i++){
    std::cout<<"Number: "<<i<<"th pose in Pose Chain"<<std::endl;
    std::cout<<"=====Ground truth Pose====="<<std::endl;
    std::cout<<poses_gt[i].matrix()<<std::endl;
    std::cout<<"=====Accumulated Pose====="<<std::endl;
    std::cout<<poses_accumulated[i].matrix()<<std::endl;
    std::cout<<"=====Optimized Pose====="<<std::endl;
    std::cout<<pose_vertices[i]->parameters().matrix()<<std::endl;
  }
  return 0;
}