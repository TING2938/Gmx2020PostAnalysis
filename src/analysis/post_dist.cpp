#include <itp/fileio>
#include <fstream>

struct Data
{
	Data(std::string fileName)
	{
		is.open(fileName);
		std::string line, buff;
		std::stringstream ss;
		std::getline(is, line);
		std::getline(is, line);
		std::getline(is, line);
		ss.str(line);
		ss >> buff >> buff >> tN1 >> buff >> buff >> tN2;
		ss.clear();
		std::getline(is, line);
		ss.str(line);
		ss >> buff >> buff >> nm1 >> buff >> buff >> nm2;
		ss.clear();
		std::getline(is, line);
		ss.str(line);
		ss >> buff >> buff >> frame;
		std::getline(is, line);
	}

	void getData(Eigen::ArrayXXi& data)
	{
		std::string line, buff;
		std::getline(is, line);
		std::stringstream ss(line);
		ss >> buff >> buff >> buff >> currTime;
		itp::loadtxt(is, data, "");
	}

	std::ifstream is;
	std::string tN1, tN2;
	int nm1, nm2, frame;
	double currTime;
};


int post_dist(int argc, char* argv[])
{
	Data data1(argv[1]), data2(argv[2]); // Nr_XX_SOL.dat, Nr_XX_Li+.dat

	/*
	** sol    0    0    1    1
	** li+    0    1    0    1
	*/
	Eigen::ArrayXd number(4);
	number.fill(0);
	int nm1 = data1.nm1;
	Eigen::ArrayXXi d1(nm1, 2), d2(nm1, 2);

	for (int i = 0; i != data1.frame; ++i)
	{
		data1.getData(d1);
		data2.getData(d2);
		ITP_ASSERT(data1.currTime == data2.currTime, "Error!");

		for (int j = 0; j != nm1; ++j)
		{
			if (d1(j, 1) == 0)
			{
				if (d2(j, 1) == 0)
					++number[0];
				else
					++number[1];
			}
			else
			{
				if (d2(j, 1) == 0)
					++number[2];
				else
					++number[3];
			}
		}
	}

	number /= data1.frame;

	std::cout << number << std::endl;

	return 0;
}

int main(int argc, char* argv[])
{
	return post_dist(argc, argv);
}
