#pragma once 


#include "solver/vertex.hh"
namespace SLAMSolver{

//Vertex Class doesn't store parameters of the Vertex
//The implementation of parameters is not exposed through Vertex
//BaseVertex stores the Parameters with ParamType
//Every vertex in solver need to derived from BaseVertex template class
//Specify:
// D: Minimal Parametrization Dimension
// ParamType: Type of Parameters
template <int D, typename ParamType>
class BaseVertex : public Vertex{
public:
  BaseVertex();
  virtual ~BaseVertex(){};
  
  ParamType parameters() const;
  
  void set_parameters(const ParamType& parameters);
  
protected:
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