#pragma once 

#include <tuple>
#include <iostream>
#include "solver/edge.hh"

namespace SLAMSolver {
template <int D, typename... VertexTypes>
class BaseEdge : public Edge{
public:
  BaseEdge();
  virtual ~BaseEdge(){};
  
  static const int NumVertices = sizeof...(VertexTypes);
  
  template <size_t Index, typename VertexType>
  void set_vertex(std::shared_ptr<VertexType> vertex_ptr);
  
  template <size_t Index>
  auto get_vertex();
    
protected:
  std::tuple<std::shared_ptr<VertexTypes>...> tuple_vertices_ptr_;
};



template <int D, typename... VertexTypes>
BaseEdge<D,VertexTypes...>::BaseEdge()
: Edge(D,NumVertices)
{
}

template <int D, typename... VertexTypes>
template <size_t Index, typename VertexType>
void BaseEdge<D,VertexTypes...>::set_vertex(std::shared_ptr<VertexType> vertex_ptr){
  //Set Template Vertex Pointer in tuple
  std::get<Index>(tuple_vertices_ptr_) = vertex_ptr;
  //Set Pointer for Vertex Interface Pointer 
  vertices_interface_ptr_[Index] = vertex_ptr;
}

template <int D, typename... VertexTypes>
template <size_t Index>
auto BaseEdge<D,VertexTypes...>::get_vertex(){
  //Get Template Vertex Pointer from tuple
  return std::get<Index>(tuple_vertices_ptr_);
}
  
  
  
  
}