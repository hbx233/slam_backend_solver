#pragma once

#include "solver/base_vertex.hh"

namespace SLAMSolver{
class InverseDepthVertex : public BaseVertex<1,double>{
  InverseDepthVertex() = default;
  void plus(const VecX &delta) override{
    parameters_ += delta(0);
  }
};
}
