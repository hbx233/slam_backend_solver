//
// Created by hbx on 2019/9/11.
//

#include "solver/sparse_schur_solver.hh"
#include "solver/base_vertex.hh"
#include "solver/base_edge.hh"
#include "solver/minimizer_levenberg_marquardt.hh"

using namespace SLAMSolver;

class DummyComplVertex : public BaseVertex<6, Eigen::Matrix<double,6,1>>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  DummyComplVertex(){}
  void plus(const VecX &delta) override{
  }
};

class DummyMargVertex : public BaseVertex<3, Eigen::Matrix<double,3,1>>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  DummyMargVertex(){}
  void plus(const VecX &delta) override{
  }
};

class DummyEdge : public BaseEdge<2,DummyComplVertex, DummyMargVertex>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  DummyEdge(){
    errors_ = Eigen::Vector2d(1,1);
    Eigen::MatrixXd jacobi0(2,6);
    jacobi0<<1,2,3,4,5,6,1,2,3,4,5,6;
    Eigen::MatrixXd jacobi1(2,3);
    jacobi1<<1,2,3,4,5,6;
    jacobians_[0] = (jacobi0);
    jacobians_[1] = (jacobi1);
    std::cout<<jacobians_[0]<<std::endl<<jacobians_[1]<<std::endl;
  }
  DummyEdge(const Eigen::Vector2d& errors, const Eigen::Matrix<double,2,6>& jacobi0, const Eigen::Matrix<double,2,3>& jacobi1)
  {
    errors_ = errors;
    jacobians_.push_back(jacobi0);
    jacobians_.push_back(jacobi1);
  }
  void compute_errors() override{
  }
  void compute_jacobians() override{
  }
};


int main(){
  std::shared_ptr<Problem> problem_ptr = std::make_shared<Problem>();
  std::vector<std::shared_ptr<DummyMargVertex>> marg_vertices;
  std::vector<std::shared_ptr<DummyComplVertex>> compl_vertices;
  //Create marginalized vertices
  for(int i = 0; i < 5; i++){
    std::shared_ptr<DummyMargVertex> marg_vertex = std::make_shared<DummyMargVertex>();
    marg_vertices.push_back(marg_vertex);
  }

  //Create complement vertices
  for(int i = 0; i < 3; i++){
    std::shared_ptr<DummyComplVertex> compl_vertex = std::make_shared<DummyComplVertex>();
    compl_vertices.push_back(compl_vertex);
  }

  //Create edges between complement vertices and marginalized vertices
  std::vector<std::shared_ptr<DummyEdge>> edges;
  std::shared_ptr<DummyEdge> edge_c0_m0 = std::make_shared<DummyEdge>();
  edge_c0_m0->set_vertex<0>(compl_vertices[0]);
  edge_c0_m0->set_vertex<1>(marg_vertices[0]);
  edges.push_back(edge_c0_m0);
  std::shared_ptr<DummyEdge> edge_c0_m1 = std::make_shared<DummyEdge>();
  edge_c0_m1->set_vertex<0>(compl_vertices[0]);
  edge_c0_m1->set_vertex<1>(marg_vertices[1]);
  edges.push_back(edge_c0_m1);
  std::shared_ptr<DummyEdge> edge_c0_m2 = std::make_shared<DummyEdge>();
  edge_c0_m2->set_vertex<0>(compl_vertices[0]);
  edge_c0_m2->set_vertex<1>(marg_vertices[2]);
  edges.push_back(edge_c0_m2);

  std::shared_ptr<DummyEdge> edge_c1_m1 = std::make_shared<DummyEdge>();
  edge_c1_m1->set_vertex<0>(compl_vertices[1]);
  edge_c1_m1->set_vertex<1>(marg_vertices[1]);
  edges.push_back(edge_c1_m1);
  std::shared_ptr<DummyEdge> edge_c1_m3 = std::make_shared<DummyEdge>();
  edge_c1_m3->set_vertex<0>(compl_vertices[1]);
  edge_c1_m3->set_vertex<1>(marg_vertices[3]);
  edges.push_back(edge_c1_m3);

  std::shared_ptr<DummyEdge> edge_c2_m2 = std::make_shared<DummyEdge>();
  edge_c2_m2->set_vertex<0>(compl_vertices[2]);
  edge_c2_m2->set_vertex<1>(marg_vertices[2]);
  edges.push_back(edge_c2_m2);
  std::shared_ptr<DummyEdge> edge_c2_m3 = std::make_shared<DummyEdge>();
  edge_c2_m3->set_vertex<0>(compl_vertices[2]);
  edge_c2_m3->set_vertex<1>(marg_vertices[3]);
  edges.push_back(edge_c2_m3);
  std::shared_ptr<DummyEdge> edge_c2_m4 = std::make_shared<DummyEdge>();
  edge_c2_m4->set_vertex<0>(compl_vertices[2]);
  edge_c2_m4->set_vertex<1>(marg_vertices[4]);
  edges.push_back(edge_c2_m4);
  //Add vertices and edges to problem
  for(auto& c_v : compl_vertices){
    problem_ptr->add_vertex(c_v);
  }

  std::set<unsigned long> marginalized_id;
  for(auto& m_v : marg_vertices){
    problem_ptr->add_vertex(m_v);
    marginalized_id.insert(m_v->id());
  }
  for(auto& e : edges){
    problem_ptr->add_edge(e);
  }

  //Create Solver
  std::shared_ptr<SparseSchurSolver> sparse_schur_solver_ptr = std::make_shared<SparseSchurSolver>(problem_ptr);
  sparse_schur_solver_ptr->set_marginalized_vertices(marginalized_id);

  //Compute Vertices start index
  std::cout<<"Compute Start Index of vertex parameters in normal equation"<<std::endl;
  sparse_schur_solver_ptr->compute_vertices_index();

  //Build Hessian and g vector
  std::cout<<"Build Hessian Matrix and g vector"<<std::endl;
  sparse_schur_solver_ptr->build_solve_structure();
  std::cout<<"Hessian: "<<std::endl;
  std::cout<<sparse_schur_solver_ptr->hessian_<<std::endl;
  std::cout<<"g: "<<std::endl;
  std::cout<<sparse_schur_solver_ptr->g_<<std::endl;

  //Compute Schur complement blocks
  std::cout<<"Compute Schur Complement Blocks"<<std::endl;
  sparse_schur_solver_ptr->compute_schur_complements();
  for(auto& schur : sparse_schur_solver_ptr->schur_complements_){
    std::cout<<"Marginalize: "<<std::endl;
    std::cout<<schur.first.start_index_<<' '<<schur.first.minimal_dim_<<std::endl;
    std::cout<<"Complement: "<<std::endl;
    for(auto& c : schur.second){
      std::cout<<c.start_index_<<' '<<c.minimal_dim_<<std::endl;
    }
  }

  //Add lambda to hessian matrix to insure non-singularity
  int size = sparse_schur_solver_ptr->hessian_.cols();
  sparse_schur_solver_ptr->hessian_ += 40 * Eigen::MatrixXd::Identity(size, size);

  std::cout<<"===delta_x from pure inverse==="<<std::endl;
  std::cout<<sparse_schur_solver_ptr->hessian_.inverse() * sparse_schur_solver_ptr->g_<<std::endl;
  //Compute delta_x_ using marginalization
  sparse_schur_solver_ptr->solve_delta_x();
  std::cout<<"===delta_x from marginalization==="<<std::endl;
  std::cout<<sparse_schur_solver_ptr->delta_x_<<std::endl;
  std::cout<<"===Solver error==="<<std::endl;
  std::cout<<sparse_schur_solver_ptr->hessian_ * sparse_schur_solver_ptr->delta_x_ - sparse_schur_solver_ptr->g_<<std::endl;
  return 0;
}
