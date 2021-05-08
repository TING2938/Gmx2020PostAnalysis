#include <vector>
#include <iomanip>
#include <itp/core>
#include <fstream>
#include <itp/timer>
#include <sstream>

int main()
{
	constexpr int N = 10000;
	fmt::print("Matrix size: {0} x {0}\n", N);
	std::vector<std::vector<double>> matrix(N);
	for (int i = 0; i < N; i++)
	{
		matrix[i].resize(N);
		for (int j = 0; j < N; j++)
		{
			matrix[i][j] = i * 1.32 + j * 4.23 + 1.42;
		}
	}

	itp::Timer timer;

	timer.start();
	std::ofstream textFile("textFile");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fmt::print(textFile, "{:15.8f} ", matrix[i][j]);
		}
		fmt::print(textFile, "\n");
	}
	textFile.close();
	timer.stop();
	fmt::print("write to text file: {} s\n", timer.span());

	timer.start();
	std::ofstream binFile("binFile", std::ios::binary);
	for (int i = 0; i < N; i++)
	{
		binFile.write((char*)(matrix.data()), sizeof(double) * N);
	}
	timer.stop();
	fmt::print("write to binary file: {} s\n", timer.span());

	std::vector<std::vector<double>> loadFromText(N);
	timer.start();
	std::string line;
	std::stringstream ss;
	std::ifstream textFileID("textFile");
	for (int i = 0; i < N; i++)
	{
		loadFromText[i].resize(N);
		std::getline(textFileID, line);
		ss.str(line);
		for (int j = 0; j < N; j++)
		{
			ss >> loadFromText[i][j];
		}
		ss.clear();
	}
	timer.stop();
	fmt::print("read from text file: {} s\n", timer.span());

	std::vector<double> loadFromBinary(N*N);
	timer.start();
	std::ifstream binFileID("binFile", std::ios::binary);
	binFileID.read((char*)(loadFromBinary.data()), sizeof(double)*N*N);

	timer.stop();
	fmt::print("read from binary file: {} s\n", timer.span());

}

