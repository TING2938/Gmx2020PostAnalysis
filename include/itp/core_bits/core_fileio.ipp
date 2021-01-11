// #include "core_fileio.hpp"
#include <fstream>

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

	template <typename T>
	Array<T, Dynamic, Dynamic> loadtxt(std::istream& is, int nrows, int ncols, std::string comments, int skiprows)
	{
		std::string line;
		T buff;
		std::stringstream ss;

		for (int i = 0; i != skiprows; ++i)
			std::getline(is, line);

		if (ncols == Dynamic)
		{
			ncols = 0;
			while (std::getline(is, line))
				if (!inner::is_comment(line, comments))
				{
					ss.str(line);
					while (ss >> buff)
						++ncols;
					break;
				}
		}
		else
		{
			std::getline(is, line);
		}

		if (nrows == Dynamic)
		{
			std::vector<T> vec;
			do
			{
				if (!inner::is_comment(line, comments))
				{
					ss.clear();
					ss.str(line);
					for (int i = 0; i != ncols; ++i)
					{
						ss >> buff;
						vec.push_back(buff);
					}
				}
			} while (std::getline(is, line));
			Array<T, Dynamic, Dynamic> res(vec.size() / ncols, ncols);
			for (int i = 0; i != res.rows(); ++i)
			{
				for (int j = 0; j != res.cols(); ++j)
				{
					res(i, j) = vec[ncols * i + j];
				}
			}
			return res;
		}

		Array<T, Dynamic, Dynamic> data(nrows, ncols);
		size_t i = 0;
		do
		{
			if (!inner::is_comment(line, comments))
			{
				ss.clear();
				ss.str(line);
				for (int j = 0; j != ncols; ++j)
					ss >> data(i, j);
				++i;
			}
		} while (i != nrows && std::getline(is, line));
		return data;
	}

	template <typename T>
	Array<T, Dynamic, Dynamic> loadtxt(std::string fileName, int nrows, int ncols, std::string comments, int skiprows)
	{
		std::ifstream is(fileName);
		ITP_ASSERT(is.is_open(), "Can not open this file: " + fileName);
		return loadtxt<T>(is, nrows, ncols, comments, skiprows);
	}

	template <typename T>
	bool loadtxt(std::istream& is, Array<T, Dynamic, Dynamic>& data, std::string comments, int skiprows)
	{
		std::string line;
		std::stringstream ss;

		for (size_t i = 0; i != skiprows; ++i)
			if (!std::getline(is, line))
			{
				return false;
			}

		for (size_t i = 0; i != data.rows(); )
		{
			if (!std::getline(is, line))
			{
				return false;
			}
			if (!inner::is_comment(line, comments))
			{
				ss.clear();
				ss.str(line);
				for (size_t j = 0; j != data.cols(); ++j)
					ss >> data(i, j);
				++i;
			}
		}
		return true;
	}

} // !namespace itp;
