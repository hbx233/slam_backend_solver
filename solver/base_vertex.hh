#pragma once 


#include "solver/vertex.hh"
namespace SLAMSolver{

/*!
 * @brief Templae base class for vertex
 * @tparam D Dimension of Minimal Parametrization
 * @tparam ParamType Type of parameter's internal representation
 */
template <int D, typename ParamType>
class BaseVertex : public Vertex{
public:
  /// @brief Constructor of BaseVertex
  BaseVertex();
  /// @brief Virtual Destructor of BaseVertex
  virtual ~BaseVertex(){};

  /// @return parameters
  ParamType parameters() const;

  /// @brief Set parameters in vertex
  void set_parameters(const ParamType& parameters);
  
protected:
  /// Storage of parameters in type of ParamType
  /// This is not exposed through Vertex base class
  ParamType parameters_;
};

template <int D, typename ParamType>
BaseVertex<D,ParamType>::BaseVertex()
: Vertex(D)
{
}

template <int D, typename ParamType>
ParamType BaseVertex<D,ParamType>::parameters() const{
  return parameters_;
}

template <int D, typename ParamType> 
void BaseVertex<D,ParamType>::set_parameters(const ParamType& parameters){
  parameters_ = parameters;
}

}