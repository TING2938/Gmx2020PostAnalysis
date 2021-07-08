#ifndef __CORE_FILEIO_HPP__
#define __CORE_FILEIO_HPP__

#include <fstream>
#include "../core"

namespace itp
{
    namespace inner
    {
        inline bool is_comment(std::string& line, std::string& comments)
        {
            auto np = line.find_first_not_of(" ");
            return np == std::string::npos || (comments.size() != 0 &&
                comments.find_first_of(line[np]) != std::string::npos);
        }
    }

    /**
     * @brief 读取文本数据文件
     * @tparam T 数据标量类型
     * @param is 文件流
     * @param nrows 行数，默认自动确定
     * @param ncols 列数，默认自动确定
     * @param comments 注释字符
     * @param skiprows 跳过开头行数
     * @return 数据矩阵
    */
    template <typename T = double>
    Eigen::Array<T, Dynamic, Dynamic> loadtxt(std::istream& is,
        int nrows = Dynamic, int ncols = Dynamic, std::string comments = "#@", int skiprows = 0)
    {
        std::string line;
        T buff;
        std::stringstream ss;

        for (int i = 0; i != skiprows; ++i)
            std::getline(is, line);

        if (ncols == Dynamic) {
            ncols = 0;
            while (std::getline(is, line))
                if (!inner::is_comment(line, comments)) {
                    ss.str(line);
                    while (ss >> buff)
                        ++ncols;
                    break;
                }
        } else {
            std::getline(is, line);
        }

        if (nrows == Dynamic) {
            std::vector<T> vec;
            do {
                if (!inner::is_comment(line, comments)) {
                    ss.clear();
                    ss.str(line);
                    for (int i = 0; i != ncols; ++i) {
                        ss >> buff;
                        vec.push_back(buff);
                    }
                }
            } while (std::getline(is, line));
            Eigen::Array<T, Dynamic, Dynamic> res(vec.size() / ncols, ncols);
            for (int i = 0; i != res.rows(); ++i) {
                for (int j = 0; j != res.cols(); ++j) {
                    res(i, j) = vec[ncols * i + j];
                }
            }
            return res;
        }

        Eigen::Array<T, Dynamic, Dynamic> data(nrows, ncols);
        int i = 0;
        do {
            if (!inner::is_comment(line, comments)) {
                ss.clear();
                ss.str(line);
                for (int j = 0; j != ncols; ++j)
                    ss >> data(i, j);
                ++i;
            }
        } while (i != nrows && std::getline(is, line));
        return data;
    }

    /**
     * @brief 读取文本数据文件
     * @tparam T 数据标量类型
     * @param fileName 文件名称
     * @param nrows 行数，默认自动确定
     * @param ncols 列数，默认自动确定
     * @param comments 注释字符
     * @param skiprows 跳过开头行数
     * @return 数据矩阵
    */
    template <typename T = double>
    Eigen::Array<T, Dynamic, Dynamic> loadtxt(std::string fileName,
        int nrows = Dynamic, int ncols = Dynamic, std::string comments = "#@", int skiprows = 0)
    {
        std::ifstream is(fileName);
        ITP_ASSERT(is.is_open(), "Can not open this file: " + fileName);
        return loadtxt<T>(is, nrows, ncols, comments, skiprows);
    }

    /**
     * @brief 读取文本数据文件并写入到data中，写入量由data的大小决定
     * @tparam T 数据标量类型
     * @param is 文件流
     * @param data 需要写入的容器
     * @param comments 注释字符
     * @param skiprows 跳过开头行数
    */
    template <typename T>
    bool loadtxt(std::istream& is, Eigen::Array<T, Dynamic, Dynamic>& data, std::string comments = "#@", int skiprows = 0)
    {
        std::string line;
        std::stringstream ss;

        for (int i = 0; i != skiprows; ++i)
            if (!std::getline(is, line)) {
                return false;
            }

        for (int i = 0; i != data.rows(); ) {
            if (!std::getline(is, line)) {
                return false;
            }
            if (!inner::is_comment(line, comments)) {
                ss.clear();
                ss.str(line);
                for (int j = 0; j != data.cols(); ++j)
                    ss >> data(i, j);
                ++i;
            }
        }
        return true;
    }

} // !namespace itp;


#endif // !__CORE_FILEIO_HPP__

/* vim: set filetype=cpp et sw=2 ts=2 ai: */
