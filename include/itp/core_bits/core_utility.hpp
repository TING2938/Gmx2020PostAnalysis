#ifndef __CORE_UTILITY_HPP__
#define __CORE_UTILITY_HPP__

namespace itp
{
    /**
     * @brief 获取本地时间日期字符串
     * @param fmt 时间日期格式，参见 https://en.cppreference.com/w/cpp/io/manip/put_time
     * @return 时间日期字符串
    */
    inline std::string localTime(const char* fmt = "%Y-%m-%d %H:%M:%S %A")
    {
        time_t t;
        char buf[500];
        t = time(NULL);
        strftime(buf, 500, fmt, localtime(&t));
        return buf;
    }

    /**
     * @brief 对数组累加
     * @tparam T 标量类型
     * @param vec 数组
     * @return 累加结果
    */
    template<typename T>
    Eigen::Array<T, -1, 1> cumsum(const Eigen::Array<T, -1, 1>& vec)
    {
        Eigen::Array<T, -1, 1> res(vec.size());
        res[0] = vec[0];
        for (int i = 1; i != vec.size(); ++i) {
            res[i] = res[i - 1] + vec[i];
        }
        return res;
    }

    /**
     * @brief 数值积分
     * @tparam S x标量类型
     * @tparam T y标量类型
     * @param x 数组x
     * @param y 数组y
     * @return 积分结果
    */
    template<typename S, typename T>
    double trapz(const Eigen::Array<S, -1, 1>& x, const Eigen::Array<T, -1, 1>& y)
    {
        ITP_ASSERT(x.size() == y.size(), "trapz: size error!");
        double res = 0.0;
        for (int i = 0; i != x.size() - 1; ++i) {
            res += (x[i + 1] - x[i]) * (y[i + 1] + y[i]);
        }
        return res / 2.0;
    }

    /**
     * @brief 累积积分
     * @tparam S x标量类型
     * @tparam T y标量类型
     * @param x  数组x
     * @param y 数组y
     * @return 累积积分数组
    */
    template<typename S, typename T>
    Eigen::ArrayXd cumtrapz(const Eigen::Array<S, -1, 1>& x, const Eigen::Array<T, -1, 1>& y)
    {
        ITP_ASSERT(x.size() == y.size(), "trapz: size error!");
        Eigen::ArrayXd res(x.size());
        res[0] = 0.0;

        for (int i = 1; i != x.size(); ++i) {
            res[i] = res[i - 1] + (x[i] - x[i - 1]) * (y[i - 1] + y[i]) / 2.0;
        }
        return res;
    }

    /**
     * @brief 数组中是否包含某元素
     * @tparam ARRAY 数组类型
     * @tparam T 元素类型
     * @param vec 数组
     * @param str 元素
     * @return 数组是否包含元素
    */
    template <typename ARRAY, typename T>
    inline bool contains(const ARRAY& vec, const T& str)
    {
        for (auto&& i : vec) {
            if (str == i)
                return true;
        }
        return false;
    }

    /**
     * @brief 在printf等函数前调用，使得字符在控制台单行滚动输出
    */
    inline void setScrollOutput()
    {
    #ifdef _WIN32
        fmt::print("\r");
    #else
        fmt::print("\r\033[k");
    #endif // _WIN32
        fflush(stdout);
    }

} // ! namespace itp

#endif // !__CORE_UTILITY_HPP__
