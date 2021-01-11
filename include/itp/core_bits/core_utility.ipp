#include "core_utility.hpp"

namespace itp
{
	inline std::string localTime(const char* fmt)
	{
		time_t t;
		char buf[500];
		t = time(NULL);
		strftime(buf, 500, fmt, localtime(&t));
		return buf;
	}

	template<typename T>
	Array<T, -1, 1> cumsum(const Array<T, -1, 1>& vec)
	{
		Array<T, -1, 1> res(vec.size());
		res[0] = vec[0];
		for (size_t i = 1; i != vec.size(); ++i)
		{
			res[i] = res[i - 1] + vec[i];
		}
		return res;
	}

	template<typename S, typename T>
	double trapz(const Array<S, -1, 1>& x, const Array<T, -1, 1>& y)
	{
		ITP_ASSERT(x.size() == y.size(), "trapz: size error!");
		double res = 0.0;
		for (size_t i = 0; i != x.size() - 1; ++i)
		{
			res += (x[i + 1] - x[i]) * (y[i + 1] + y[i]);
		}
		return res / 2.0;
	}

	template<typename S, typename T>
	ArrayXd cumtrapz(const Array<S, -1, 1>& x, const Array<T, -1, 1>& y)
	{
		ITP_ASSERT(x.size() == y.size(), "trapz: size error!");
		ArrayXd res(x.size());
		res[0] = 0.0;

		for (size_t i = 1; i != x.size(); ++i)
		{
			res[i] = res[i - 1] + (x[i] - x[i - 1]) * (y[i - 1] + y[i]) / 2.0;
		}
		return res;
	}

	template <typename ARRAY, typename T>
	inline bool contains(const ARRAY& vec, const T& str)
	{
		for (auto&& i : vec)
		{
			if (str == i)
				return true;
		}
		return false;
	}

	inline void setScrollOutput()
	{
#ifdef _WIN32
		fmt::print("\r");
#else
		fmt::print("\r\033[k");
#endif // _WIN32
		fflush(stdout);
	}
}
