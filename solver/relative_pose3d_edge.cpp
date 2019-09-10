#include "solver/relative_pose3d_edge.h"

namespace SLAMSolver{
RelativePose3DEdge::RelativePose3DEdge(const Sophus::SE3 &pose0_SE3_pose1_meas)
:pose0_SE3_pose1_meas_(pose0_SE3_pose1_meas)
{
}

RelativePose3DEdge::
void RelativePose3DEdge::compute_errors() {
  //Fetch pose from vertex
  Sophus::SE3 frame_SE3_pose0 = get_vertex<0>()->parameters();
  Sophus::SE3 frame_SE3_pose1 = get_vertex<1>()->parameters();
  //Compute Error between measurement and estimation
  Sophus::SE3 error_SE3 =  pose0_SE3_pose1_meas_.inverse() * frame_SE3_pose0.inverse() * frame_SE3_pose1;
  //Separately use lie algebra of rotation and translation instead of
  //using full lie algebra of transformation
  errors_ = error_SE3.log();
}

void RelativePose3DEdge::compute_jacobians() {
  //Fetch pose from vertex
  Sophus::SE3 frame_SE3_pose0 = get_vertex<0>()->parameters();
  Sophus::SE3 frame_SE3_pose1 = get_vertex<1>()->parameters();
  //Compute jacobians for pose0

}
}