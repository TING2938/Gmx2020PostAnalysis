#include <itp/core>
#include <thread>
#include <vector>

void func()
{
	fmt::print("Test subprogram\n");
}

class TestClass
{
public:
	void classFunc()
	{
		fmt::print("class func test");
	}
};

int main(int argc, char* argv[])
{
	
#define RED(str) "\033[31m" str "\033[0m"
#define YEL(str) "\033[33m" str "\033[0m"

	itp::Getopt getopt(argc, argv, "Test for Getopt class");

	fmt::print(RED("aaa")YEL("bbb"));

	int a1 = 0, a2 = 3;
	double b1 = 3.5, b2 = 5.6;
	string s1{ "s1" }, s2{ "s2" };
	vecd v1{ 2, 4, 5 };

	getopt(a1, "-a1", false, "a1");
	getopt(b1, "-b1", false, "b1");
	getopt(s1, "-s1", false, "s1");
	getopt.getArray(v1, "-v1", false, "v1"); 
	getopt.finish();

}