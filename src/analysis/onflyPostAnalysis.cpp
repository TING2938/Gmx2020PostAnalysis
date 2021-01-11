#include <itp/core>
#include <itp/getopt>
#include <itp/fileio>

void getShape(std::ifstream& file, const std::string& sign, int& nrow, int& ncol)
{
	std::string line;
	std::stringstream ss;
	double tmp;

	std::getline(file, line);
	std::getline(file, line);
	ss.str(line);
	ncol = 0;
	while (ss >> tmp) 
	{
		ncol++;
	}
	nrow = 1;
	while (std::getline(file, line))
	{
		if (line == sign)
		{
			break;
		}
		nrow++;
	}
	file.clear();
	file.seekg(std::ios::beg);
}

int main(int argc, char** argv)
{
	std::string fnm = "densMNC.dat";
	std::string sign = "20080513\tgfeng\t";
	std::string output = "newdensMNC.dat";
	int freq = 1000;
	double dt = 0.002;
	double beginTime = 0;
	double endTime = 0;

	itp::Getopt getopt(argc, argv);
	getopt(freq, "-freq", true, "freq, userint2 in .mdp file");
	getopt(dt, "-dt", true, "dt in .mdp file");
	getopt(beginTime, "-b", false, "begin time (ps)");
	getopt(endTime, "-e", false, "end time (ps)");
	getopt(fnm, "-f", false, "input file");
	getopt(sign, "-s", false, "sign");
	getopt(output, "-o", false, "output file name");
	getopt.finish();

	fmt::print("freq: {}\n", freq);
	fmt::print("dt: {}\n", dt);
	fmt::print("beginTime: {}\n", beginTime);
	fmt::print("endTine: {}\n", endTime);
	fmt::print("sign: {}\n", sign);
	fmt::print("input file name: {}\n", fnm);
	fmt::print("output file name: {}\n", output);
	
	dt = dt * freq;
	size_t beginStep = beginTime / dt;
	size_t endStep = endTime / dt;
	if (endStep == 0)
	{
		endStep = -1;
	}

	int nrow = 0;
	int ncol = 0;
	auto file = std::ifstream(fnm);
	getShape(file, sign, nrow, ncol);
	
	fmt::print("ncol: {}\n", ncol);
	fmt::print("nbin: {}\n", nrow);

	Eigen::ArrayXXd data(nrow, ncol), allData(nrow, ncol);
	allData.fill(0);
	data.fill(0);
	std::string line;

	for (size_t i = 0; i < beginStep * (nrow + 2); i++)
	{
		if (!std::getline(file, line))
		{
			fmt::print("Error, no more frame!\n");
			std::exit(-1);
		}
	}

	fmt::print("Analysis ...\n");
	size_t totFrame = 0;
	for (size_t i = beginStep; i <= endStep; i++)
	{
		if (!itp::loadtxt(file, data))
		{
			break;
		}
		allData += data;
		std::getline(file, line);
		std::getline(file, line);
		totFrame++;
	}
	fmt::print("number of frame: {}\n", totFrame);
	fmt::print("Analysis done!\n");

	allData /= totFrame;
	
	std::ofstream outputfile(output);
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < ncol; j++)
		{
			fmt::print(outputfile, "{:15.8e}\t", allData(i, j));
		}
		fmt::print(outputfile, "\n");
	}
}

