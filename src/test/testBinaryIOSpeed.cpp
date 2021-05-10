#include <vector>
#include <iomanip>
#include <itp/core>
#include <fstream>
#include <itp/timer>
#include <sstream>

using real = float;

int main()
{
	constexpr int N = 5000;
	fmt::print("Matrix size: {0} x {0}\n", N);
	std::vector<std::vector<real>> matrix(N);
	for (int i = 0; i < N; i++)
	{
		matrix[i].resize(N);
		for (int j = 0; j < N; j++)
		{
			matrix[i][j] = i * i * 1.32 + j * 4.23 + 1.42;
		}
	}

	itp::Timer timer;


	std::ofstream binFile("binFile", std::ios::binary);
	timer.start();
	for (int i = 0; i < N; i++)
	{
		binFile.write((char*)(matrix[i].data()), sizeof(real) * N);
	}
	timer.stop();
	binFile.close();
	fmt::print("write to binary file: {} s\n", timer.span());

	std::vector<real> loadFromBinary(N*N);
	std::ifstream binFileID("binFile", std::ios::binary);
	timer.start();
	binFileID.read((char*)(loadFromBinary.data()), sizeof(real)*N*N);
	timer.stop();
	binFileID.close();
	fmt::print("read from binary file: {} s\n", timer.span());
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (loadFromBinary[i * N + j] != matrix[i][j])
			{
				fmt::print("Not equal: [{}][{}]: {} and {}\n", i, j, loadFromBinary[i * N + j], matrix[i][j]);
				std::exit(-1);
			}
		}
	}
	fmt::print("All equal\n");
}

