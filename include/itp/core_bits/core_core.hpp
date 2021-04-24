#ifndef __CORE_CORE_H__
#define __CORE_CORE_H__

#define EIGEN_NO_DEBUG

#include <string>
#include <vector>
#include <array>
#include <iostream>

#include <Eigen/Dense>

#include <fmt/format.h>
#include <fmt/ostream.h>
using namespace fmt::literals;


// assert
#ifdef ITP_NDEBUG
#define ITP_ASSERT(expression, comments) ((void)0)
#else
#define ITP_ASSERT(expression, comments) (void)(              \
			(!!(expression)) ||                                   \
			((std::cerr << __FILE__ << "(" << __LINE__ << ")\n" << comments << std::endl, exit(-2), 0))   \
	)
#endif


// math
namespace itp
{
	constexpr long double pi = 3.14159265358979323846264338327950288L;
	constexpr long double exp1 = 2.71828182845904523536028747135266250L;
}

// type define
namespace itp
{
	using Eigen::Dynamic;

	template <typename T>
	using Mat = Eigen::Array<T, Dynamic, Dynamic>;

	template <typename T>
	using Vec = Eigen::Array<T, Dynamic, 1>;

	using boxd = Mat<Eigen::Array3d>;

	using matd = Eigen::ArrayXXd;
	using mati = Eigen::ArrayXXi;

	using vecd = Eigen::ArrayXd;
	using veci = Eigen::ArrayXi;
}

#endif // !__CORE_CORE_H__

/* vim: set filetype=cpp et sw=2 ts=2 ai: */
