#include <itp/core>

#include <format>

#define Print(...) std::cout << std::format(__VA_ARGS__)

int main()
{
	Eigen::Vector3d vec(1, 3, 5);
	std::cout << vec << std::endl;

	std::cout << std::format("a is {}\n", 42);

	Print("b is {} {}\n", 43, 42);
	
}

