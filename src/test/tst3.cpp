#include <itp/getopt>
#include <fmt/color.h>
#include <vector>
#include <string>
#include <itp/logger>
using namespace Eigen;

int main(int argc, char** argv)
{
	Eigen::ArrayXXd arr(12, 1);
	arr << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
	
	Eigen::ArrayXXd B(12, 3);
	B << arr, arr, arr;


	fmt::print("{}\n", B);
	

	


}
