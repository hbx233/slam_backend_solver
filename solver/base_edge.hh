#pragma once 

#include <tuple>
#include <iostream>
#include "solver/edge.hh"

namespace SLAMSolver {
/*!
 * @brief The variadic template base class for all error edges
 * @tparam D Dimension of errors
 * @tparam VertexTypes Types of Vertices the edge connected to
 */
template <int D, typename... VertexTypes>
class BaseEdge : public Edge{
public:
  /// @brief Constructor, initialize the Edge base class
  /// with errors dimension and number of vertices
  BaseEdge();

  /// @brief Virtual destructor
  virtual ~BaseEdge(){};

  /// @brief Number of vertices specified in the template list
  static const int NumVertices = sizeof...(VertexTypes);

  /*!
   * @brief Set vertex pointer in the tuple
   * @tparam Index Index of the vertex want to set
   * @tparam VertexType Type of the vertex, this can be deducted, no need to specify
   * @param vertex_ptr The vertex pointer want to set
   */
  template <size_t Index, typename VertexType>
  void set_vertex(std::shared_ptr<VertexType> vertex_ptr);

  /*!
   * @brief Get vertex pointer in the tuple
   * @tparam Index Index of the vertex want to get
   * @return vertex pointer with deducted type
   */
  template <size_t Index>
  auto get_vertex();
    
protected:
  /// A tuple containing all vertices' pointer that the edge connected to
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