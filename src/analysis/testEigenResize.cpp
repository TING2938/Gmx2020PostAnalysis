#include <itp/core>
#include <itp/fileio>

int main()
{
    auto data = itp::loadtxt("test.txt");
    Eigen::ArrayXd tmp = data.row(0).middleCols(3, 6);
    tmp.resize(3, 2);
    std::cout << "data: \n" << data 
        << "\ntmp: \n" << tmp << std::endl;
    std::cout << tmp.rowwise().sum() << std::endl;
}