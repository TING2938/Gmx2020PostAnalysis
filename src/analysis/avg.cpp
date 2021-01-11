#include <itp/getopt>	
#include <itp/fileio>

int avg(int argc, char* argv[])
{
	itp::Getopt getopt(argc, argv, "Calculate statistical information of data.");
	std::string fileName;
	std::vector<int> col{ 2 };
	int skiprows = 0;
	constexpr double inf = std::numeric_limits<double>::infinity();
	double begTime = inf;
	double endTime = inf;

	getopt.getFixPos(fileName, 1, true, "file name to analysis");
	getopt.getArray(col, "-c", false, "Column to analysis.");
	std::sort(col.begin(), col.end());
	getopt(begTime, "-b", false, "Begin position.");
	getopt(endTime, "-e", false, "End position.");
	getopt(skiprows, "-skip", false, "Skipped line.");
	getopt.finish();

	auto data = itp::loadtxt(argv[1], Eigen::Dynamic, Eigen::Dynamic, "#@", skiprows);

	auto ll = data.rows();
	double t1 = data(0, 0);
	double t2 = data(ll - 1, 0);

	if (begTime == inf)
		begTime = t1;
	if (endTime == inf)
		endTime = t2;

	if (argc >= 3 && isdigit(argv[2][0]))
		begTime = atof(argv[2]);
	if (argc >= 4 && isdigit(argv[2][0]) && isdigit(argv[3][0]))
		endTime = atof(argv[3]);

	size_t b = static_cast<size_t>((begTime - t1) * (ll - 1) / (t2 - t1));
	size_t e = static_cast<size_t>((endTime - t1) * (ll - 1) / (t2 - t1)) + 1;

	auto ncol = col.size();
	Eigen::ArrayX4d result(ncol + 1, 4); // mean, std, max, min;

	Eigen::ArrayXd sum(e - b);
	sum.fill(0);


	Eigen::ArrayXd temp;
	for (size_t i = 0; i != ncol; ++i)
	{
		temp = data.middleRows(b, e - b).col(col[i] - 1);
		result(i, 0) = temp.mean();
		result(i, 1) = std::sqrt((temp - temp.mean()).square().sum() / (temp.size() - 1));
		result(i, 2) = temp.maxCoeff();
		result(i, 3) = temp.minCoeff();
		sum += temp;
	}

	result(ncol, 0) = sum.mean();
	result(ncol, 1) = std::sqrt((temp - temp.mean()).square().sum() / (sum.size() - 1));
	result(ncol, 2) = sum.maxCoeff();
	result(ncol, 3) = sum.minCoeff();

	printf("\nStatistics from %.3f to %.3f, over %zd points\n", data(b, 0), data(e - 1, 0), e - b);
	printf("\n%8s%13s%13s%13s%13s\n", "Column", "Average", "Std", "Max", "Min");
	fmt::print("{}\n", std::string(13 * 4 + 8, '-'));

	for (int i = 0; i != ncol; ++i)
		printf("%8d%13.3f%13.3f%13.3f%13.3f\n", col[i], result(i, 0), result(i, 1), result(i, 2), result(i, 3));
	if (ncol != 1)
	{
		fmt::print("{}\n", std::string(13 * 4 + 8, '-'));
		printf("%8s%13.3f%13.3f%13.3f%13.3f\n", "Totol", result(ncol, 0), result(ncol, 1), result(ncol, 2), result(ncol, 3));
	}
	printf("\n");
	return 0;
}

int main(int argc, char* argv[])
{
	avg(argc, argv);
}
