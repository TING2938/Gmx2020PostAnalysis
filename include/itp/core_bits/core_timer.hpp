#ifndef _CORE_TIMER_H_
#define _CORE_TIMER_H_

#include <chrono>
#include <iostream>
#include <string>

/*
usage: 
itp::Timer timer;
timer.start();
// dosomething();
timer.stop();
fmt::print("It takes {} seconds\n", timer.span());
*/

namespace itp
{
	/**
	 * @brief 计时器类
	*/
	class Timer
	{
	public:
		using nanoseconds = std::chrono::duration<double, std::nano>;   // ns
		using microseconds = std::chrono::duration<double, std::micro>; // us
		using milliseconds = std::chrono::duration<double, std::milli>; // ms
		using seconds = std::chrono::duration<double>;                  // s
		using minutes = std::chrono::duration<double, std::ratio<60>>;  // min
		using hours = std::chrono::duration<double, std::ratio<3600>>;  // h
		
	private:
		std::chrono::steady_clock::time_point _begin, _end;
		
	public:
		Timer() : _begin(std::chrono::steady_clock::time_point()),
			_end(std::chrono::steady_clock::time_point())
		{
		}

		/**
		 * @brief 开始计时
		*/
		void start()
		{
			_begin = std::chrono::steady_clock::now();
		}

		/**
		 * @brief 结束计时
		*/
		void stop()
		{
			_end = std::chrono::steady_clock::now();
		}

		/**
		 * @brief 计时跨度
		 * @return 从开始到结束所用时间（默认：秒）
		*/
		template <typename Duration = seconds>
		double span()
		{
			return std::chrono::duration_cast<Duration>(_end - _begin).count();
		}

	};

} // namespace

#endif // _CORE_TIMER_H_
