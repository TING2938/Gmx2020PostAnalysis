#include <iostream>
#include <itp/core>
#include <fstream>
#include <vector>
#include <string>
#include <itp/timer>
using namespace std;

int main(int argc, char** argv)
{

	// for read
	int len = stoi(argv[2]);

	itp::Timeit timeit;

	std::vector<float> output(100);

	timeit.start();
	for (int i = 0; i < 100; i++) {
		ifstream ifile(argv[1], ios::binary);
		vector<float> vec(len * len);
		ifile.read((char*)(vec.data()), len * len * sizeof(float));
		ifile.close();
		output[i] = vec[0];
	}
	timeit.stop();

	fmt::print("{}: spend {} ms\n", output[0], timeit.span<itp::Timeit::milliseconds>());
}
