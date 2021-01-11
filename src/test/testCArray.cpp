#include <itp/core>

void func(int (&vec)[3], int vec2[3])
{
	fmt::print("{}\n", sizeof(vec));
	vec[2] = 3;
	vec2[2] = 100;
}

int main()
{
	int vec[3] = {2, 4, 5};
	int vec2[3] = {1, 2, 3};
	func(vec, vec2);
}