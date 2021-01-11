#include <filesystem>
#include <itp/core>

int main()
{
	printf("%ws\n", std::filesystem::current_path().c_str());
}
