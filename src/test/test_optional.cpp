#include <optional>
#include <iostream>

std::optional<int> func(bool a)
{
	if (a)
	{
		return 5;
	}
	else
	{
		return std::nullopt;
	}
}

int main()
{
	if (auto ret = func(true); ret.has_value())
	{
		std::cout << ret.value() << std::endl;
	}
	else
	{
		std::cout << "no value." << std::endl;
	}
}