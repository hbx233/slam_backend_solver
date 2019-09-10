#pragma once 

#include "solver/eigen_types.h"

namespace SLAMSolver{
class Vertex {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  Vertex(const int minimal_dimension);

  virtual ~Vertex();

  int minimal_dimension() const{return minimal_dimension_;}

  unsigned long id() const { return id_; }

  virtual void plus(const VecX &delta) = 0;

  void set_fixed(bool fixed = true) {
      fixed_ = fixed;
  }

  bool is_fixed() const { return fixed_; }

protected:
  const int minimal_dimension_;
  unsigned long id_;

  bool fixed_ = false;   //True if the Vertex is fixed
};
}