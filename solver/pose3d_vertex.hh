#pragma once
#include "solver/base_vertex.hh"
#include "sophus/se3.hpp"

namespace SLAMSolver {
class Pose3DVertex : public BaseVertex<6,Sophus::SE3>{
public:
  Pose3DVertex(){};
  //Override the update function for 3D pose
  void plus(const VecX &delta) override;
};
}
