#include <itp/core>

#include <format>
#include <map>
#include <vector>

#define Print(...) std::cout << std::format(__VA_ARGS__)

std::vector<std::string> splitString(const std::string& str)
{
    std::vector<std::string>          result;
    std::string::const_iterator       currPos = str.begin();
    const std::string::const_iterator end = str.end();
    while (currPos != end) {
        while (currPos != end && std::isspace(*currPos)) {
            ++currPos;
        }
        const std::string::const_iterator startPos = currPos;
        while (currPos != end && !std::isspace(*currPos)) {
            ++currPos;
        }
        if (startPos != end) {
            result.emplace_back(startPos, currPos);
        }
    }
    return result;
}


int main()
{
	Eigen::Vector3d vec(1, 3, 5);
	std::cout << vec << std::endl;

	std::cout << std::format("a is {}\n", 42);

	Print("b is {} {}\n", 43, 42);
    std::map<std::string, std::vector<std::string>> map;
    auto ret = splitString("Value 32 65");
    map[ret[0]] = std::vector<std::string>(ret.begin()+1, ret.end());

    std::vector<int> va = { 1, 2, 3 };
    va.clear();
    va.push_back(2);
    va.clear();
    va.push_back(3);
    va.push_back(6);
	
}

