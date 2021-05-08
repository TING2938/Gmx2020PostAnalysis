#include <Eigen/Dense>
#include <fmt/format.h>
#include <itp/timer>
#include <iostream>
#include <cstdio>
#include <algorithm>

extern"C" {
	//LU decomoposition of a general matrix
	void dgetrf_(int* M, int* N, double* A, int* lda, int* IPIV, int* INFO);
	//generate inverse of a matrix given its LU decomposition
	void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

void inverse_Lapack(double* A, double* B, int N)
{
	std::copy_n(A, N*N, B);
	int* IPIV = new int[N + 1];
	int LWORK = N * N;
	double* WORK = new double[LWORK];
	int INFO;
	dgetrf_(&N, &N, B, &N, IPIV, &INFO);
	dgetri_(&N, B, &N, IPIV, WORK, &LWORK, &INFO);
	delete[] IPIV;
	delete[] WORK;
}

int main()
{
	int N;
	fmt::print("Please input number: ");
	std::cin >> N;

	Eigen::MatrixXd mat = Eigen::MatrixXd::Random(N, N) * 10;
	Eigen::MatrixXd inv = mat;
	itp::Timer timer;

	timer.start();
	inv = mat.inverse();
	timer.stop();
	fmt::print("{}, Eigen time : {}s\n", inv(N / 2, N / 2), timer.span());

	timer.start();
	inverse_Lapack(mat.data(), inv.data(), N);
	timer.stop();
	fmt::print("{}, Lapack time: {}s\n", inv(N / 2, N / 2), timer.span());
}