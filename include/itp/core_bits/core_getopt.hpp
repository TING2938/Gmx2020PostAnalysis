#ifndef __CORE_GETOPT_H__
#define __CORE_GETOPT_H__

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>

namespace itp
{
	/**
	 * @brief 获取命令行选项
	*/
	class Getopt
	{
	public:
		
		Getopt() = default;
		
		/**
		 * @brief 构造函数
		 * @param gc argc
		 * @param gv argv
		 * @param msg 程序说明
		 * @return
		*/
		Getopt(int gc, char** gv, std::string msg = "") : argc(gc), _mainMessage(msg)
		{
			for (int i = 0; i != gc; i++)
			{
				argv.push_back(gv[i]);
			}
			auto p = findArgs("-h");
			_printHelpMain = (p == argv.begin() + 1);
			_printHelpSub = (argc == 3 && p == argv.begin() + 2);

		}
		
		Getopt(std::string str, std::string msg="") : _mainMessage(msg)
		{
			std::stringstream ss(str);
			std::string tmp;
			while (ss >> tmp)
			{
				argv.push_back(tmp);
			}
			argc = argv.size();
			auto p = findArgs("-h");
			_printHelpMain = (p == argv.begin() + 1);
			_printHelpSub = (argc == 3 && p == argv.begin() + 2);
		}
		
		void setCommand(std::string str, std::string msg="")
        {
            _mainMessage = msg;
            std::stringstream ss(str);
            std::string tmp;
            while (ss >> tmp)
            {
                argv.push_back(tmp);
            }
            argc = argv.size();
            auto p = findArgs("-h");
            _printHelpMain = (p == argv.begin() + 1);
            _printHelpSub = (argc == 3 && p == argv.begin() + 2);
        }

        void setCommand(int gc, char** gv, std::string msg = "")
        {
            argc = gc;
            _mainMessage = msg;
            for (int i = 0; i != gc; i++)
            {
                argv.push_back(gv[i]);
            }
            auto p = findArgs("-h");
            _printHelpMain = (p == argv.begin() + 1);
            _printHelpSub = (argc == 3 && p == argv.begin() + 2);
        }
		
		/**
		 * @brief 获取标量类型命令行参数
		 * @tparam T 标量类型
		 * @param x 输入数据
		 * @param str 命令行参数指示字符串
		 * @param required 是否是必需参数
		 * @param msg 参数说明
		*/
		template<typename T>
		void operator()(T& x, std::string str, bool required = false, std::string msg = "")
		{
			if (addHelpInfo(required, str, x, msg))
				return;

			auto p = findArgs(str);
			if (p != argv.end())
			{
				std::stringstream ss{ *(p + 1) };
				ss >> x;
			}
			else checkRequired(required, str);
		}

		/**
		 * @brief 获取布尔类型命令行参数
		 * @param x 输入数据
		 * @param str 命令行参数指示字符串
		 * @param required 是否是必需参数
		 * @param msg 参数说明
		*/
		void operator()(bool& x, std::string str, bool required = false, std::string msg = "")
		{
			if (addHelpInfo(required, str, x, msg))
				return;

			auto p = findArgs(str);
			if (p != argv.end())
				x = true;
			else checkRequired(required, str);
		}

		/**
		 * @brief 获取数组类型命令行参数
		 * @tparam T 标量类型
		 * @param x 输入数据
		 * @param str 命令行参数指示字符串
		 * @param required 是否是必需参数
		 * @param msg 参数说明
		*/
		template <typename T>
		void getArray(std::vector<T>& x, std::string str, bool required = false, std::string msg = "")
		{
			std::stringstream ss;
			ss << "[";
			for (auto ii : x)
			{
				ss << ii << ", ";
			}
			if (!x.empty())
				ss << "\b\b]";
			else
				ss << "]";
			if (addHelpInfo(required, str, ss.str(), msg))
				return;

			auto b = findArgs(str);
			if (b != argv.end())
			{
				x.clear();
				auto e = std::find_if(b + 1, argv.end(),
					[](const std::string& chr) { return chr[0] == '-'; });
				T tmp;
				for (size_t i = 1; i != e - b; ++i)
				{
					std::stringstream ss{ *(b + i) };
					ss >> tmp;
					x.push_back(tmp);
				}
			}
			else checkRequired(required, str);
		}

		/**
		 * @brief 获取固定位置命令行参数，从1开始表示第一个命令行
		 * @tparam T 标量类型
		 * @param x 输入数据
		 * @param pos 参数位置
		 * @param required 是否是必需参数
		 * @param msg 参数说明
		*/
		template <typename T>
		void getFixPos(T& x, int pos, bool required = false, std::string msg = "")
		{
			if (addHelpInfo(required, "pos: " + std::to_string(pos), x, msg))
				return;

			std::stringstream ss{ argv[pos] };
			ss >> x;
		}

		/**
		 * @brief 添加子程序
		 * @tparam Func 子程序执行的函数类型
		 * @tparam ...Valty 子程序执行函数的输入参数类型
		 * @param str 命令行参数指示字符串
		 * @param msg 子程序说明
		 * @param func 子程序执行的函数
		 * @param ...val 子程序执行函数的输入参数
		*/
		template <typename Func, typename... Valty>
		void addSubProgram(std::string str, std::string msg, Func func, Valty&&... val)
		{
			if (_printHelpMain)
			{
				_subprogramInfo += ("     " + str + "     " + msg + "\n");
				return;
			}

			if (str == argv[1])
			{
				func(std::forward<Valty>(val)...);
			}
		}

		/**
		 * @brief 获取完所有参数后调用
		*/
		void finish()
		{
			if (_printHelpMain)
			{
				if (!_mainMessage.empty())
				{
					printf("%s\n", _mainMessage.c_str());
				}
				printf("Command line option:\n");
				printf("%s\n\n", toString().c_str());

				if (!(_requiredInfo.empty() && _optionalInfo.empty()))
				{
					printf("\033[31m%11s%14s%36s\033[0m\n", "Option", "Value", "Description");
					printf("%s\n", std::string(61, '-').c_str());

					if (!_requiredInfo.empty())
					{
						printf("\033[33m(Required)\033[0m\n");
						printf("%s", _requiredInfo.c_str());
					}
					if (!_optionalInfo.empty())
					{
						printf("\033[33m(Optional)\033[0m\n");
						printf("%s", _optionalInfo.c_str());
					}
					printf("%s\n\n", std::string(61, '-').c_str());
				}

				if (!_subprogramInfo.empty())
				{
					printf("Sub Program Message:\n");
					printf("\033[31m%15s     %s\033[0m\n", "Function", "Description");
					printf("%s\n", std::string(56, '-').c_str());

					printf("%s", _subprogramInfo.c_str());
					printf("%s\n\n", std::string(56, '-').c_str());
				}
			}
			if (_printHelpSub)
			{
				if (!_mainMessage.empty())
				{
					printf("%s\n", _mainMessage.c_str());
				}
				printf("Command line option:\n");
				printf("%s\n\n", toString().c_str());

				if (!(_requiredInfo.empty() && _optionalInfo.empty()))
				{
					printf("\033[31m%15s%25s     %s\033[0m\n", "Option", "Value", "Description");
					printf("%s\n", std::string(56, '-').c_str());

					if (!_requiredInfo.empty())
					{
						printf("\033[33m(Required)\033[0m\n");
						printf("%s", _requiredInfo.c_str());
					}
					if (!_optionalInfo.empty())
					{
						printf("\033[33m(Optional)\033[0m\n");
						printf("%s", _optionalInfo.c_str());
					}
					printf("%s\n\n", std::string(56, '-').c_str());
				}
			}
			if (_printHelpMain || _printHelpSub)
			{
				std::exit(0);
			}
		}

		/**
		 * @brief 将命令行参数转化为字符串
		 * @return 命令行参数字符串
		*/
		std::string toString()
		{
			std::string str(argv[0]);
			for (size_t i = 1; i != argc; ++i)
			{
				str += " ";
				str += argv[i];
			}
			return str;
		}

	private:
		void checkRequired(bool required, std::string str)
		{
			if (required)
			{
				printf("Command Line Error: ");
				printf("\"%s\"", str.c_str());
				printf(" option NOT Found, type \"-h\" for help.\nexit!\n");
				std::exit(-2);
			}
		}

		template <typename T>
		bool addHelpInfo(bool required, std::string str, const T& x, std::string msg)
		{
			if (_printHelpMain || _printHelpSub)
			{
				std::stringstream ss; 
				ss << std::setiosflags(std::ios::left) << std::string(5, ' ')  << std::setw(15) << str << std::setw(30) << x << msg << "\n";
				if (required)
				{
					_requiredInfo += ss.str();
				}
				else
				{
					_optionalInfo += ss.str();
				}
				return true;
			}
			return false;
		}

		/*! \brief
		 *  return argument position, or (argv+argc) if not found.
		 */
		std::vector<std::string>::iterator findArgs(std::string str)
		{
			return std::find(argv.begin()+1, argv.end(), str);
		}

	private:
		int    argc;
		std::vector<std::string> argv;
		bool _printHelpMain;
		bool _printHelpSub;
		std::string _mainMessage;
		std::string _requiredInfo;
		std::string _optionalInfo;
		std::string _subprogramInfo;  // for sub program;

	}; // class Getopt;

}


#endif // !__CORE_GETOPT_H__
