#ifndef __ITP_CORE_LOGGER_HPP__
#define __ITP_CORE_LOGGER_HPP__

#include <fstream>
#include <string>
#include <itp/utility>

/*
usage: 
int i = 4;
double j = 2.3;
std::string k = "this is a string";

itp::Logger log("log_onfly.log");
LOG(log, i);
LOG(log, j);
LOG(log, k);
*/

namespace itp
{
	class Logger
	{
	public:
		Logger(std::string fnm)
		{

			_file.open(fnm, std::ios::app);
			_file << "\nLog file opened at: " << localTime() << std::endl;
		}

		~Logger()
		{
			_file.close();
		}

		template <typename T>
		void print(const T& var, std::string name, int line, std::string path)
		{
			_file << "[" << path << ":" << line << "] " << name << ": " << var << std::endl;
		}

#define LOG(log, v) log.print(v, #v, __LINE__, __FILE__)

	private:
		std::ofstream _file;
	};
}

#endif // !__ITP_CORE_LOGGER_HPP__
