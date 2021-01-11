#ifndef __CORE_CONFIG_H__
#define __CORE_CONFIG_H__

// #define EIGEN_NO_DEBUG

// #define ITP_NDEBUG

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

#endif // !__CORE_CONFIG_H__

/* vim: set filetype=cpp et sw=2 ts=2 ai: */
