#include <iostream>
#include <itp/core>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

int main()
{
	// for write
	vector<double> vec;
	for (int i = 0; i < 10000; i++)
	{
		vec.push_back(i * 0.5);
	}
	int len = vec.size();
	string str = "length if vector is: ";

	ofstream ofile("binaryFile.dat", ios::binary);
	ofile.write((char*)str.data(), str.size() * sizeof(char));
	ofile.write((char*)&len, sizeof(len));
	ofile.write((char*)(vec.data()), len * sizeof(double));
	ofile.close();

	ofstream textFile("textFile.dat");
	for (int i = 0; i < 100; i++)
	{
		for (int j = 0; j < 100; j++)
		{
			fmt::print(textFile, "{:12.8e} ", vec[100 * i + j]);
		}
		fmt::print(textFile, "\n");
	}
	textFile.close();

	// for read
	/*
	ifstream ifile("binaryFile.dat", ios::binary);
	int len;
	ifile.read((char*)&len, sizeof(int));
	vector<double> vec(len);
	ifile.read((char*)(vec.data()), len * sizeof(double));
	for (auto&& i : vec)
	{
		cout << i << endl;
	}
	*/
}