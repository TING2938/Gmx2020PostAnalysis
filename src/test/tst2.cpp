#define CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>

#include <itp/core>
#include <vector>

int func(int argc, char* argv[])
{
	veci x = {1, 0, 3, 4, 7};
	veci y = {5, 0, 6, 4, 0};

	itp::cumsum(x).print("test for cumsum x:");
	itp::cumsum(y).print("test for cumsum y:");

	x.setUnion(y).print();
	x.setIntersection(y).print();

	mati mat(2, 2, {5, 8, 0, 6});
	mati mat1(2, 2, {5, 2, 0, 6});
	itp::find(mat == mat1).print();
	auto b = itp::find(mat > 5);
	b.print();
	itp::cumtrapz(x, y).print();
	x.flip().print();

	veci v0{1, 3, 2, 5};
	v0.cumprod().print();

	auto dif = v0.diff();
	dif.print();
	veci v1 = {1, 3, 4, 3, 2, 2, 10, 5, 8, 4, 9, 8, 6, 3, 5, 6, 7, 7, 7, 8, 9, 9, 9, 9, 10};
	veci v2 = itp::arange(5);
	veci v3 = v1.setDifference(v0);
	v3.print();
	auto aa = v1.contains(9);

	return 0;
}

int main(int argc, char* argv[])
{
	func(argc, argv);

	_CrtDumpMemoryLeaks();
}

