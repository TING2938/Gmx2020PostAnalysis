#include <Eigen/Sparse>
#include <iostream>

int main()
{
	Eigen::SparseVector<double> vec(10);
	vec.insert(2) = 23;
	vec.insert(5) = 13.5;
	std::cout << vec << std::endl;

	Eigen::SparseMatrix<double> mat(10, 10);
	mat.insert(2, 4) = 2.4;
	mat.insert(5, 1) = 3.53;
	std::cout << mat << std::endl;

}