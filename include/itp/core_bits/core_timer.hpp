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
		typedef std::chrono::steady_clock::time_point   TimePoint;
		typedef std::chrono::duration<double>           Duration;
		typedef std::chrono::steady_clock               SteadyClock;
	private:
		TimePoint _begin;
		TimePoint _end;
	public:
		Timer() : _begin(TimePoint()), _end(TimePoint())
		{
		}

		/**
		 * @brief 开始计时
		*/
		void start()
		{
			_begin = SteadyClock::now();
		}

		/**
		 * @brief 结束计时
		*/
		void stop()
		{
			_end = SteadyClock::now();
		}

		/**
		 * @brief 计时跨度
		 * @return 从开始到结束所用时间（秒）
		*/
		double span()
		{
			Duration span = std::chrono::duration_cast<Duration>(_end - _begin);
			return span.count();
		}

	};

} // namespace

#endif // _CORE_TIMER_H_
