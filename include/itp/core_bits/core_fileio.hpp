#ifndef __CORE_FILEIO_HPP__
#define __CORE_FILEIO_HPP__


namespace itp
{
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
	Array<T, Dynamic, Dynamic> loadtxt(std::istream& is,
		int nrows = Dynamic, int ncols = Dynamic, std::string comments = "#@", int skiprows = 0);

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
	Array<T, Dynamic, Dynamic> loadtxt(std::string fileName,
		int nrows = Dynamic, int ncols = Dynamic, std::string comments = "#@", int skiprows = 0);

	/**
	 * @brief 读取文本数据文件并写入到data中，写入量由data的大小决定
	 * @tparam T 数据标量类型
	 * @param is 文件流
	 * @param data 需要写入的容器
	 * @param comments 注释字符
	 * @param skiprows 跳过开头行数
	*/
	template <typename T>
	bool loadtxt(std::istream& is, Eigen::Array<T, Dynamic, Dynamic>& data, std::string comments = "#@", int skiprows = 0);

} // !namespace itp;

#include "core_fileio.ipp"

#endif // !__CORE_FILEIO_HPP__

/* vim: set filetype=cpp et sw=2 ts=2 ai: */
