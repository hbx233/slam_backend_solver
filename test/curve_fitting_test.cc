#include <iostream>
#include <random>
#include "solver/problem.hh"
#include "solver/base_vertex.hh"
#include "solver/base_edge.hh"

using namespace SLAMSolver;
using namespace std;


//Curve Fitting Vertex
//Minimal Parametrization Dimension is 3
//Use Eigen::Vector3d to represent internal Parameter
class CurveFittingVertex: public SLAMSolver::BaseVertex<3,Eigen::Vector3d>
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  CurveFittingVertex(){}

  void plus(const VecX &delta) override;
};

//Override the
void CurveFittingVertex::plus(const VecX &delta)
{
    parameters_(0) += delta(0);
    parameters_(1) += delta(1);
    parameters_(2) += delta(2);
}

// Edge for Curve Fitting
// One Edge Just Connect to one CurveFittingVertex
class CurveFittingEdge: public BaseEdge<1,CurveFittingVertex>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    CurveFittingEdge( double x, double y ):x_(x),y_(y) {}

    virtual void compute_errors() override
    {
        //First get current estimated curve parameters
        //No need to perform casting
        Eigen::Vector3d abc = get_vertex<0>()->parameters();
        //Compute error_ vector using current estimated parameters
        errors_(0) = std::exp( abc(0)*x_*x_ + abc(1)*x_ + abc(2) ) - y_;
    }
    virtual void compute_jacobians() override
    {
        //Get current estimated curve parameters
        Eigen::Vector3d abc = get_vertex<0>()->parameters();
        double exp_y = std::exp( abc(0)*x_*x_ + abc(1)*x_ + abc(2) );
        //Compute jacobian matrix
        Eigen::Matrix<double, 1, 3> jaco_abc;
        jaco_abc << x_ * x_ * exp_y, x_ * exp_y , 1 * exp_y;
        jacobians_[0] = jaco_abc;
    }
public:
    //Measurent of x and y
    double x_;
    double y_;
};

int main()
{
    double a=2.0, b=3.0, c=1.0;//Ground Truth
    int N = 200;
    double w_sigma= 1.;

    std::default_random_engine generator;
    std::normal_distribution<double> noise(0.,w_sigma);

    // 构建 problem
    Problem problem;
    shared_ptr< CurveFittingVertex > vertex(new CurveFittingVertex());
    vertex->set_parameters(Eigen::Vector3d (0,0,0));
    problem.add_vertex(vertex);

    for (int i = 0; i < N; ++i) {
        double x = i/100.;
        double n = noise(generator);
        double y = std::exp( a*x*x + b*x + c ) + n;
        //Create Edge(measurement) and add it to problem 
        shared_ptr< CurveFittingEdge > edge(new CurveFittingEdge(x,y));
        edge->set_vertex<0>(vertex);

        problem.add_edge(edge);
    }

    std::cout<<"\nTest CurveFitting start..."<<std::endl;
    problem.solve(30);
    std::cout << vertex->parameters() << std::endl;
    return 0;
}


