#include <itp/core>

int main()
{
	Eigen::Matrix3d m1;
	m1 << 1.223, 1.335, 65.4,
		24.78, 83.22234, 23.222,
		21.12, 432.3, 2322.34;

	fmt::print("{}", m1);
}