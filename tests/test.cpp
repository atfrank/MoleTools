#include <Eigen/Dense>
#include <iostream>
#include <cstdlib>

int main (){
	Eigen::MatrixXd cmp; //Nx3 Mobile Coordinate Matrix
  Eigen::MatrixXd ref; //Nx3 Stationary Coordinate Matrix
  Eigen::Matrix3d R; //3x3 Correlation Matrix
  //Eigen::JacobiSVD<Eigen::MatrixXd> svd(R, Eigen::ComputeThinU | Eigen::ComputeThinV); //SVD Matrix
  Eigen::MatrixXd U; //Equivalent to Bosco's V
  Eigen::MatrixXd S;
  Eigen::MatrixXd V; //Equivalent to Bosco's Wt
  double UdVd; //U determinant * V determinant
  unsigned int i;

  cmp.resize(20,3);
  ref.resize(20,3);

  for (i=0; i< 20; i++){
    cmp.row(i) << i, -i, i+4;
		ref.row(i) << i*i, i*i+2, -(i*i);
  }

  R=ref.transpose()*cmp;

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(R, Eigen::ComputeThinU | Eigen::ComputeThinV); //SVD Matrix

  U=svd.matrixU();
  S=svd.singularValues();
  V=svd.matrixV();

  UdVd=U.determinant()*V.determinant();

  //if (UdVd > 0.0){
  //  std::cout << "HERE\n";
  //}

  //std::cout << U << std::endl;

  std::cout << U*V << std::endl;
}
