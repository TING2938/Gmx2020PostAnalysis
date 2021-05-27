// https://github.com/99x/timercpp

#ifndef _CORE_TIMER_H_
#define _CORE_TIMER_H_

#include <chrono>
#include <iostream>
#include <string>
#include <atomic>
#include <thread>

namespace itp
{
	class Timeit
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
		Timeit() : _begin(std::chrono::steady_clock::time_point()),
			_end(std::chrono::steady_clock::time_point())
		{}


		void start()
		{
			_begin = std::chrono::steady_clock::now();
		}

		// 结束计时
		void stop()
		{
			_end = std::chrono::steady_clock::now();
		}

		// @brief 计时跨度
		// @return 从开始到结束所用时间（默认：秒）
		template <typename Duration = seconds>
		double span()
		{
			return std::chrono::duration_cast<Duration>(_end - _begin).count();
		}

		// @brief 打印出计时跨度（默认：秒）
		template <typename Duration = seconds>
		void printSpan(const std::string& front = "", const std::string& back = "")
		{
			std::cout
				<< front
				<< std::chrono::duration_cast<Duration>(_end - _begin).count()
				<< back;
		}

	};

	class Timer
	{
		std::atomic<bool> active{ true };

	public:

		/**
		 * @brief
		 * @tparam _Fn
		 * @tparam ..._Arg
		 * @param delay ms
		 * @param func function
		*/
		template <typename _Fn>
		void setTimeout(int delay, _Fn&& func)
		{
			active = true;
			std::thread t([=]() {
				if (!active.load()) return;
				std::this_thread::sleep_for(std::chrono::milliseconds(delay));
				if (!active.load()) return;
				func();
				});
			t.detach();
		}

		/**
		 * @brief
		 * @tparam _Fn
		 * @tparam ..._Arg
		 * @param interval ms
		 * @param func  function
		*/
		template <typename _Fn>
		void setInterval(int interval, _Fn&& func)
		{
			active = true;
			std::thread t([=]() {
				while (active.load()) {
					std::this_thread::sleep_for(std::chrono::milliseconds(interval));
					if (!active.load()) return;
					func();
				}
				});
			t.detach();
		}

		void stop()
		{
			active = false;
		}
	};

} // namespace

#endif // _CORE_TIMER_H_
