#include <vector>
#include <iomanip>
#include <itp/core>
#include <fstream>
#include <itp/timer>
#include <sstream>

int main()
{
	std::vector<std::vector<double>> matrix(1000);
	for (int i = 0; i < 1000; i++)
	{
		matrix[i].resize(1000);
		for (int j = 0; j < 1000; j++)
		{
			matrix[i][j] = i * 1.32 + j * 4.23 + 1.42;
		}
	}

	itp::Timer timer;

	timer.start();
	std::ofstream textFile("textFile");
	for (int i = 0; i < 1000; i++)
	{
		for (int j = 0; j < 1000; j++)
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
	for (int i = 0; i < 1000; i++)
	{
		binFile.write((char*)(matrix.data()), sizeof(double) * 1000);
	}
	timer.stop();
	fmt::print("write to binary file: {} s\n", timer.span());

	std::vector<std::vector<double>> loadFromText(1000);
	timer.start();
	std::string line;
	std::stringstream ss;
	std::ifstream textFileID("textFile");
	for (int i = 0; i < 1000; i++)
	{
		loadFromText[i].resize(1000);
		std::getline(textFileID, line);
		ss.str(line);
		for (int j = 0; j < 1000; j++)
		{
			ss >> loadFromText[i][j];
		}
		ss.clear();
	}
	timer.stop();
	fmt::print("read from text file: {} s\n", timer.span());

	std::vector<double> loadFromBinary(1000*1000);
	timer.start();
	std::ifstream binFileID("binFile", std::ios::binary);
	binFileID.read((char*)(loadFromBinary.data()), sizeof(double)*1000*1000);

	timer.stop();
	fmt::print("read from binary file: {} s\n", timer.span());

}

