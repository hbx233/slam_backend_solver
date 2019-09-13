//
// Created by hbx on 2019/9/11.
//
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
#include "solver/pose3d_vertex.hh"
#include "solver/reprojection_edge.hh"
#include "solver/sparse_schur_solver.hh"
#include "solver/problem.hh"
#include "solver/minimizer_levenberg_marquardt.hh"

using namespace SLAMSolver;

//Test for monocular sparse bundle adjustment
struct Frame {
  Frame(const Sophus::SE3& world_SE3_camera) : world_SE3_camera_(world_SE3_camera){}
  Sophus::SE3 world_SE3_camera_;
  std::unordered_map<int, Eigen::Vector3d> keypoint_id_map_keypoint_;
};


void get_simulation_data(std::vector<Frame>& frames, std::unordered_map<int, Eigen::Vector3d>& points){
  int feature_nums = 100;
  int frame_nums = 5;

  double radius = 8;

  for (int n = 0; n < frame_nums; ++n) {
    double theta = n * 2 * M_PI / (frame_nums * 4);

    Sophus::SO3 R(0,0,theta);
    Eigen::Vector3d t = Eigen::Vector3d(radius * cos(theta) - radius, radius * sin(theta), 1 * sin(2 * theta));
    frames.push_back(Frame({R, t}));
  }

  std::default_random_engine generator;
  std::normal_distribution<double> noise_pdf(0., 1. / 1000.);  // 2pixel / focal
  for (int j = 0; j < feature_nums; ++j) {
    //Generate one feature point
    std::uniform_real_distribution<double> xy_rand(-4, 4.0);
    std::uniform_real_distribution<double> z_rand(4., 8.);

    Eigen::Vector3d point_world(xy_rand(generator), xy_rand(generator), z_rand(generator));
    points[j] = point_world;

    //generate observation - projected point with noise
    for (int i = 0; i < frame_nums; ++i) {
      Eigen::Vector3d point_camera = frames[i].world_SE3_camera_.inverse() * point_world;
      point_camera /= point_camera[2];
      //perturb observation with noise
      point_camera[0] += noise_pdf(generator);
      point_camera[1] += noise_pdf(generator);
      frames[i].keypoint_id_map_keypoint_.insert(std::make_pair(j, point_camera));
    }
  }
}

int main() {
  std::vector <Frame> frames;
  std::unordered_map<int, Eigen::Vector3d> points;
  get_simulation_data(frames, points);

  std::normal_distribution<double> noise_pdf(0., 1);
  std::default_random_engine generator;

  //create sparse bundle adjustment problem
  std::shared_ptr<Problem> problem_ptr = std::make_shared<Problem>();
  //Create pose and keypoint vertices
  //First provide with ground truth state
  std::unordered_map<int, std::shared_ptr<Point3DVertex>> point_id_map_point_vertex;
  std::vector<std::shared_ptr<Pose3DVertex>> pose_vertices;
  //cache the vertex value before optimization
  std::unordered_map<int, Eigen::Vector3d> point_id_map_point_init;
  std::vector<Sophus::SE3> poses_init;

  std::set<unsigned long> marginalized_id_set;
  for(auto& point : points){
    std::shared_ptr<Point3DVertex> point_ptr = std::make_shared<Point3DVertex>();
    Eigen::Vector3d point_value = point.second;
    //perturb the point
    point_value(0) += noise_pdf(generator);
    point_value(1) += noise_pdf(generator);
    point_value(2) += noise_pdf(generator);

    point_ptr->set_parameters(point_value); //set to ground truth point position in world
    point_id_map_point_init[point.first] = point_value;
    marginalized_id_set.insert(point_ptr->id()); //add landmark point id to marginalized set
    point_id_map_point_vertex[point.first] = point_ptr; //store the vertex pointer
    //Add vertices to problem
    problem_ptr->add_vertex(point_ptr);
  }

  //Create pose vertices
  for(int i = 0; i < frames.size(); i++){
    std::shared_ptr<Pose3DVertex> pose_ptr = std::make_shared<Pose3DVertex>();
    if(i < 2){
      //monocular, fix the first two frames to set a scale
      //also ensure the optimizer won't drift in nullspace
      pose_ptr->set_fixed();
      pose_ptr->set_parameters(frames[i].world_SE3_camera_);
    } else{
      //perturb the frame after the first two frames
      pose_ptr->set_parameters(frames[1].world_SE3_camera_);
    }
    //pose_ptr->set_fixed();
    //pose_ptr->set_parameters(frames[0].world_SE3_camera_);
    pose_vertices.push_back(pose_ptr);
    poses_init.push_back(pose_ptr->parameters());
    //Add vertices to problem
    problem_ptr->add_vertex(pose_ptr);
  }

  //Create edges between 3D keypoint and pose
  for(int i = 0; i < frames.size(); i++){
    //for every observation
    for(auto& obs : frames[i].keypoint_id_map_keypoint_){
      std::shared_ptr<ReprojectionEdge> edge_ptr = std::make_shared<ReprojectionEdge>(obs.second.segment(0,2));
      //Set keypoint vertex and pose vertex in ReprojectionEdge
      edge_ptr->set_vertex<0>(point_id_map_point_vertex[obs.first]);
      edge_ptr->set_vertex<1>(pose_vertices[i]);
      //Add edge to problem
      problem_ptr->add_edge(edge_ptr);
    }
  }

  //Create SparseSchurSolver
  std::shared_ptr<SparseSchurSolver> solver_ptr = std::make_shared<SparseSchurSolver>(problem_ptr);
  //add marginalized id set to solver for marginalization sparse schur operation
  solver_ptr->set_marginalized_vertices(marginalized_id_set);

  //Create LM minimizer
  LevenbergMarguardtConfig lm_config;
  MinimizerLevenbergMarquardt lm_minimizer(solver_ptr, lm_config);
  lm_minimizer.minimize(50);

  //Output the solver result
  std::cout<<"=====Camera Pose Optimization Result====="<<std::endl;
  for(int i = 0; i < frames.size(); i++){
    std::cout<<"Camera "<<i<<std::endl;
    std::cout<<"Ground Truth: "<<std::endl;
    std::cout<<frames[i].world_SE3_camera_.matrix()<<std::endl;
    std::cout<<"Before Optimization: "<<std::endl;
    std::cout<<poses_init[i].matrix()<<std::endl;
    std::cout<<"Optimized: "<<std::endl;
    std::cout<<pose_vertices[i]->parameters().matrix()<<std::endl;
  }

  std::cout<<"=====Keypoint 3D position Optimization Result====="<<std::endl;
  for(auto& point : points){
    std::cout<<"Keypoint "<<point.first<<std::endl;
    std::cout<<"Ground Truth: "<<std::endl;
    std::cout<<point.second<<std::endl;
    std::cout<<"Before Optimization: "<<std::endl;
    std::cout<<point_id_map_point_init[point.first]<<std::endl;
    std::cout<<"Optimized: "<<std::endl;
    std::cout<<point_id_map_point_vertex[point.first]->parameters()<<std::endl;
  }
}
