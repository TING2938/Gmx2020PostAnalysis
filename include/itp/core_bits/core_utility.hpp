#ifndef __CORE_UTILITY_HPP__
#define __CORE_UTILITY_HPP__

namespace itp
{
	/**
	 * @brief 获取本地时间日期字符串
	 * @param fmt 时间日期格式，参见 https://en.cppreference.com/w/cpp/io/manip/put_time
	 * @return 时间日期字符串
	*/
	inline std::string localTime(const char* fmt = "%Y-%m-%d %H:%M:%S %A");

	/**
	 * @brief 对数组累加
	 * @tparam T 标量类型
	 * @param vec 数组
	 * @return 累加结果
	*/
	template<typename T>
	Array<T, -1, 1> cumsum(const Array<T, -1, 1>& vec);

	/**
	 * @brief 数值积分
	 * @tparam S x标量类型
	 * @tparam T y标量类型
	 * @param x 数组x
	 * @param y 数组y
	 * @return 积分结果
	*/
	template<typename S, typename T>
	double trapz(const Array<S, -1, 1>& x, const Array<T, -1, 1>& y);

	/**
	 * @brief 累积积分
	 * @tparam S x标量类型
	 * @tparam T y标量类型
	 * @param x  数组x
	 * @param y 数组y
	 * @return 累积积分数组
	*/
	template<typename S, typename T>
	ArrayXd cumtrapz(const Array<S, -1, 1>& x, const Array<T, -1, 1>& y);

	/**
	 * @brief 数组中是否包含某元素
	 * @tparam ARRAY 数组类型
	 * @tparam T 元素类型
	 * @param vec 数组
	 * @param str 元素
	 * @return 数组是否包含元素
	*/
	template <typename ARRAY, typename T>
	inline bool contains(const ARRAY& vec, const T& str);

	/* set  scroll output, only valid for the next line.*/

	/**
	 * @brief 在printf等函数前调用，使得字符在控制台单行滚动输出
	*/
	inline void setScrollOutput();

}

#include "core_utility.ipp"

#endif // !__CORE_UTILITY_HPP__
