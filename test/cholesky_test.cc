//
// Created by hbx on 2019/9/15.
//
#include "Eigen/Eigen"
#include <iostream>
using namespace Eigen;
using namespace std;
int main(){

  MatrixXd A(6,6);
  A << 9.00018, 9.44428e-05, 5.78167e-05, -0.00372673, 2.27808, 28.7904,
       9.44428e-05,9.00012 ,-7.85597e-05  ,   -2.27889 , -0.00207573    , -47.9074,
  5.78167e-05, -7.85597e-05   ,    9.0002   ,   -28.793    ,   47.912 , -0.00126306,
  -0.00372673  ,   -2.27889 ,     -28.793   ,   20207.2   ,  -341.178   ,   35.2242,
  2.27808 , -0.00207573   ,    47.912  ,   -341.178  ,    20576.2   ,    20.679,
  28.7904 ,    -47.9074 , -0.00126306  ,    35.2242   ,    20.679    ,  20777.7;
  cout << "The matrix A is" << endl << A << endl;
  LLT<MatrixXd> lltOfA(A); // compute the Cholesky decomposition of A
  MatrixXd L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
  MatrixXd U = lltOfA.matrixU();
  // The previous two lines can also be written as "L = A.llt().matrixL()"
  cout << "The Cholesky factor L is" << endl << L << endl;
  cout << "The Cholesky factor U is" << endl << U << endl;
  cout << "To check this, let us compute L * L.transpose()" << endl;
  cout << L * L.transpose() << endl;
  cout << " U * U.transpose() " << endl;
  cout << U.transpose() * U<< endl;
  cout << "This should equal the matrix A" << endl;
  return 0;

}
