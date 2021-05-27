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

	itp::Timer timer;

	std::vector<float> output(100);

	timer.start();
	for (int i = 0; i < 100; i++) {
		ifstream ifile(argv[1], ios::binary);
		vector<float> vec(len * len);
		ifile.read((char*)(vec.data()), len * len * sizeof(float));
		ifile.close();
		output[i] = vec[0];
	}
	timer.stop();

	fmt::print("{}: spend {} ms\n", output[0], timer.span<itp::Timer::milliseconds>());
}