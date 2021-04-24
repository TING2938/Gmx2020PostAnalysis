#include <itp/core>
#include <itp/timer>
#include <Eigen/Dense>
#include <array>

double** creatMatrix(int m, int n)
{
	double** ret = new double* [m];
	for (int i = 0; i < m; i++)
	{
		ret[i] = new double[n];
	}
	return ret;
}

std::vector<std::vector<double>> creatStdMatrix(int m, int n)
{
	std::vector<std::vector<double>> ret(m);
	for (auto&& r : ret)
		r.resize(n);
	return ret;
}

int main()
{
	itp::Timer time;
	int N = 100000000;
	{
		time.start();
		auto mat1 = creatStdMatrix(N, 3);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				mat1[i][j] = 3.25 * i + 22.64 * j;
			}
		}

		time.stop();
		fmt::print("vector value: {}, time: {} ms\n", mat1[N/2][1], time.span<itp::Timer::milliseconds>());
	}

	{
		time.start();
		auto mat3 = creatMatrix(N, 3);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				mat3[i][j] = 3.25 * i + 22.64 * j;
			}
		}

		time.stop();
		fmt::print("new value: {}, time: {} s\n", mat3[N / 2][1], time.span());
	}

	{
		time.start();
		std::vector<std::array<double, 3>> mat2(N);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				mat2[i][j] = 3.25 * i + 22.64 * j;
			}
		}

		time.stop();
		fmt::print("std vector array value: {}, time: {} ms\n", mat2[N / 2][1], time.span<itp::Timer::milliseconds>());
	}
	
	{
		using rvec = double[3];
		time.start();
		auto mat4 = new rvec[N];

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				mat4[i][j] = 3.25 * i + 22.64 * j;
			}
		}

		time.stop();
		fmt::print("new plain value: {}, time: {} s\n", mat4[N / 2][1], time.span());
	}

	{
		time.start();
		Eigen::ArrayXXd mat5(N, 3);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				mat5(i, j) = 3.25 * i + 22.64 * j;
			}
		}

		time.stop();
		fmt::print("Eigen value: {}, time: {} s\n", mat5(N / 2, 1), time.span());
	}

}
