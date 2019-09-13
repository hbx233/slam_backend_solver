#pragma once
#include "solver/base_edge.hh"
#include "pose3d_vertex.hh"
#include "point3d_vertex.hh"
namespace SLAMSolver {
//Represent the relative pose measurement from pose at 1 to pose at 0
//pose0_SE3_pose1
class RelativePose3DEdge : public BaseEdge<2,Point3DVertex,Pose3DVertex>{
public:
  RelativePose3DEdge(const Sophus::SE3& pose0_SE3_pose1_meas);
  Eigen::Matrix<double,6,6> compute_jacobian_of_SE3();
  void compute_errors() override;
  void compute_jacobians() override;
private:
  Sophus::SE3 pose0_SE3_pose1_meas_;
};
}