#include <itp/getopt>
#include <fmt/color.h>
#include <array>
#include <vector>
#include <string>
#include <map>
#include <itp/logger>
#include <itp/timer>
#include <numeric>
#include <functional>

int nbin2 = 3000000;

double efunc(double r)
{
	constexpr double ewaldcoeff = 1.255521;
	constexpr double one4pieps0 = 138.935458;
	constexpr double eta = 19.79;
	return std::erfc(ewaldcoeff * r) / r - std::erfc(eta * r) / r;
}

class InterpFuncTable
{
	double low, up;
	std::vector<double> dbin, limit, startx, allTable;
	std::vector<int> nbin;

public:
	void setLowUp(double Low, double Up)
	{
		low = Low;
		up = Up;
	}
	void setDbin(const std::vector<double>& Dbin)
	{
		dbin = Dbin;
	}
	void setLimit(const std::vector<double>& Limit)
	{
		limit = Limit;
	}
	void makeTable(std::function<double(double)> func, bool bCheck = false)
	{
		double LOW = low*0.95;
		double UP = up * 1.05;
		nbin.resize(dbin.size() + 1, 0);
		startx.resize(dbin.size(), 0);
		if (limit.size() != dbin.size())
		{
			printf("Error! limit.size() != dbin.size()\n");
			std::exit(-1);
		}
		double tmpl = LOW, tmpu = 0;
		for (int i = 0; i < dbin.size(); i++)
		{
			startx[i] = tmpl;
			if (tmpl > limit[i])
				continue;
			tmpu = std::min(UP, limit[i]);
			nbin[i + 1] = nbin[i] + int((tmpu - tmpl) / dbin[i]) + 50;
			for (int j = 0; j < int((tmpu - tmpl) / dbin[i]) + 50; j++)
			{
				allTable.push_back(func(tmpl + dbin[i] * j));
			}
			tmpl = tmpu;
		}

		printf("Nbin: ");
		for (auto i : nbin)
		{
			printf("%d ", i);
		}
		printf("\n");

		// for check
		if (bCheck)
		{
			int nbin2 = 300000;
			double dbin2 = (up - low) / nbin2;
			std::vector<double> y1(nbin2), y2(nbin2), dy(nbin2);

			for (int i = 0; i < nbin2; i++)
			{
				y1[i] = func(low + dbin2 * i);
				y2[i] = this->at(low + dbin2 * i);
				dy[i] = std::abs(y1[i] - y2[i]);
			}
			auto argmax = std::max_element(dy.begin(), dy.end()) - dy.begin();
			fmt::print("max dy[{}]: {}\n", argmax, dy[argmax]);
			auto dd = this->at(low + dbin2 * argmax);
			printf("sum of dy : %f\n", std::accumulate(dy.begin(), dy.end(), 0.0));
		}
	}

	double at(double r)
	{
		for (int i = startx.size() - 1; i >= 0; i--)
		{
			if (r >= startx[i])
			{
				int me = (r - startx[i]) / dbin[i];
				return allTable[me + nbin[i]] + (r - (startx[i] + dbin[i] * me)) * 
					(allTable[me + nbin[i] + 1] - allTable[me + nbin[i]]) / dbin[i];
			}
		}
		return 0.0;
	}
};

class InterpFuncTable_log
{
	double low, up, logUp, logLow;
	int nbin;
	double dlog;
	std::vector<double> allX, allTable;

public:
	void setParam(double Low, double Up, double Dlog)
	{
		logUp = std::log(Up);
		logLow = std::log(Low) + Dlog;
		low = Low;
		up = std::exp(logUp);
		dlog = Dlog;
		nbin = (logUp - logLow) / dlog;
		fmt::print("nbin: {}\n", nbin);
	}

	void makeTable(std::function<double(double)> func, bool bCheck = false)
	{
		allX.resize(nbin, 0);
		allTable.resize(nbin, 0);
		for (int i = 0; i < nbin; i++)
		{
			allX[i] = std::exp(logLow + i * dlog);
			allTable[i] = func(allX[i]);
		}

		// for check
		if (bCheck)
		{
			int nbin2 = 3000000;
			double dbin2 = (up / 1.05 - low / 0.95) / nbin2;
			std::vector<double> y1(nbin2), y2(nbin2), dy(nbin2);

			itp::Timeit timeit;

			timeit.start();
			for (int i = 0; i < nbin2; i++)
			{
				y1[i] = func(low + dbin2 * i);
			}
			timeit.stop();
			timeit.printSpan("func: ", " s\n");

			timeit.start();
			for (int i = 0; i < nbin2; i++)
			{
				y2[i] = this->at(low + dbin2 * i);
			}
			timeit.stop();
			timeit.printSpan("at: ", " s\n");

			for (int i = 0; i < nbin2; i++)
			{
				dy[i] = std::abs(y1[i] - y2[i]);
			}

			auto argmax = std::max_element(dy.begin(), dy.end()) - dy.begin();
			fmt::print("max dy[{}]: {}\n", argmax, dy[argmax]);
			auto dd = this->at(low + dbin2 * argmax);
			printf("%f, sum of dy : %f\n", dd, std::accumulate(dy.begin(), dy.end(), 0.0));
		}
	}

	double at(double r)
	{
		int me = (std::log(r) - logLow) / dlog;
		return allTable[me] + (r - allX[me]) * (allTable[me + 1] - allTable[me]) / (allX[me + 1] - allX[me]);
	}
};

void getLimitLiner()
{
	fmt::print("#-------liner\n");
	double LOW = 0.003;
	double UP = 3;

	int nbin = 30000;
	double low = LOW * 0.95;
	double up = UP * 1.05;
	double dbin = (up - low) / nbin;

	std::vector<double> allTable(nbin);
	for (int i = 0; i < nbin; i++)
	{
		allTable[i] = efunc(low + dbin * i);
	}
	fmt::print("{} {}\n", low+dbin*(nbin), low + dbin*(nbin-1));

	auto findTable = [&](double r) {
		int me = (r - low) / dbin;
		return allTable[me] + (r - low - dbin * me) * (allTable[me + 1] - allTable[me]) / dbin;
	};
	
	/* ===================== */
	double dbin2 = (UP - LOW) / nbin2;
	std::vector<double> y1(nbin2), y2(nbin2), dy(nbin2);

	itp::Timeit timeit;

	timeit.start();
	for (int i = 0; i < nbin2; i++)
	{
		y1[i] = efunc(LOW + dbin2 * i);
	}
	timeit.stop();
	timeit.printSpan("efunc: ", "s\n");

	timeit.start();
	for (int i = 0; i < nbin2; i++)
	{
		y2[i] = findTable(LOW + dbin2 * i);
	}
	timeit.stop();
	timeit.printSpan("table: ", "s\n");

	for (int i = 0; i < nbin2; i++)
	{
		dy[i] = std::abs(y1[i] - y2[i]);
	}
	double maxdy = *std::max_element(dy.begin(), dy.end());
	fmt::print("max dy: {}\n", maxdy);
	fmt::print("{}\n", std::accumulate(dy.begin(), dy.end(), 0.0));
}

void getLimitLog()
{
	fmt::print("#-------log\n");
	double LOW = 0.003;
	double UP = 3;

	int nbin = 10000;
	double low = LOW * 0.95;
	double up = UP * 1.05;
	double logUp = std::log(up);
	double logLow = std::log(low);
	double dlog = (logUp - logLow) / nbin;

	std::vector<double> allX(nbin), allTable(nbin);
	for (int i = 0; i < nbin; i++)
	{
		allX[i] = std::exp(std::log(low) + i*dlog);
		allTable[i] = efunc(allX[i]);
	}

	auto findTable = [&](double r) {
		int me = (std::log(r) - logLow) / dlog;
		return allTable[me] + (r - allX[me]) * (allTable[me + 1] - allTable[me]) / (allX[me+1]-allX[me]);
	};

	double dbin2 = (UP - LOW) / nbin2;
	std::vector<double> y1(nbin2), y2(nbin2), dy(nbin2);

	itp::Timeit timeit;

	timeit.start();
	for (int i = 0; i < nbin2; i++)
	{
		y1[i] = efunc(LOW + dbin2 * i);
	}
	timeit.stop();
	timeit.printSpan("efunc: ", "s\n");

	timeit.start();
	for (int i = 0; i < nbin2; i++)
	{
		y2[i] = findTable(LOW + dbin2 * i);
	}
	timeit.stop();
	timeit.printSpan("table: ", "s\n");

	for (int i = 0; i < nbin2; i++)
	{
		dy[i] = std::abs(y1[i] - y2[i]);
	}
	double maxdy = *std::max_element(dy.begin(), dy.end());
	fmt::print("max dy: {}\n", maxdy);
	fmt::print("{}\n", std::accumulate(dy.begin(), dy.end(), 0.0));
}

void getLimitCut()
{
	fmt::print("#-------Cut\n");
	double LOW = 0.003;
	double UP = 1.2;

	double low = LOW * 0.95;
	double up = UP * 1.05;

	std::vector<double> dbin = { 2e-5, 1e-4, 0.00025 };
	std::vector<double> limit = {0.15, 1   , 10000.0 };

	std::vector<int> nbin(dbin.size()+1, 0);
	std::vector<double> startx(dbin.size(), 0);
	std::vector<double> allTable;
	double tmpl = low, tmpu = 0;
	for (int i = 0; i < dbin.size(); i++)
	{
		startx[i] = tmpl;
		if (tmpl > limit[i])
			continue;
		if (UP < limit[i])
		{
			tmpu = up;
		}
		else {
			tmpu = limit[i];
		}
		nbin[i+1] = nbin[i] + int((tmpu - tmpl) / dbin[i])+100;
		for (int j = 0; j < int((tmpu - tmpl) / dbin[i])+100; j++)
		{
			allTable.push_back(efunc(tmpl + dbin[i] * j));
		}
		tmpl = tmpu;
	}

	auto findTable = [&](double r) {
		for (int i = startx.size()-1; i >= 0; i--)
		{
			if (r > startx[i])
			{
				int me = (r - startx[i]) / dbin[i];
				if (me + nbin[i] + 1 > allTable.size())
				{
					fmt::print("i: {}\n", i);
				}
				return allTable[me+nbin[i]] + (r - (startx[i] + dbin[i] * me)) * (allTable[me + nbin[i] + 1] - allTable[me + nbin[i]]) / dbin[i];
			}
		}
		return 0.0;
	};

	/* ===================== */

	double dbin2 = (UP - LOW) / nbin2;
	std::vector<double> y1(nbin2), y2(nbin2), dy(nbin2);

	itp::Timeit timeit;

	timeit.start();
	for (int i = 0; i < nbin2; i++)
	{
		y1[i] = efunc(LOW + dbin2 * i);
	}
	timeit.stop();
	timeit.printSpan("efunc: ", "s\n");

	timeit.start();
	for (int i = 0; i < nbin2; i++)
	{
		y2[i] = findTable(LOW + dbin2 * i);
	}
	timeit.stop();
	timeit.printSpan("table: ", "s\n");

	for (int i = 0; i < nbin2; i++)
	{
		dy[i] = std::abs(y1[i] - y2[i]);
	}
	auto maxdy = std::max_element(dy.begin(), dy.end()) - dy.begin();

	fmt::print("max dy[{}]: {}\n", maxdy, dy[maxdy]);
	auto dd = findTable(LOW + dbin2 * maxdy);
	auto dd1 = efunc(LOW + dbin2 * maxdy);
	fmt::print("{} \n{}\n", dd - dd1, std::accumulate(dy.begin(), dy.end(), 0.0));


}

int main()
{
	fmt::print("\n\nfor class errorFunc\n");
	InterpFuncTable_log errorFunc;
	errorFunc.setParam(0.003 * 0.003, 3 * 3, 3e-4);
	float epsilon_r = 2;
	float eta = 19.79;
	float ewaldcoeff = 1.255521;
	float alpha = ewaldcoeff;
	#define ONE_4PI_EPS0 138.935458

	errorFunc.makeTable([epsilon_r, alpha, eta](double r2) {
		double sqrtM_PI = std::sqrt(itp::pi);
		double alpha2 = alpha * alpha;
		double eta = 19.79;
		double eta2 = eta * eta;
		double r = std::sqrt(r2);
		double r3 = r * r2;
		double f_r = 2 * alpha * std::exp(-alpha2 * r2) / (sqrtM_PI * r2)
			+ erfc(alpha * r) / (r3)
			-(2 * eta * std::exp(-eta2 * r2) / (sqrtM_PI * r2) + erfc(eta * r) / (r3));
		return epsilon_r * f_r * ONE_4PI_EPS0;
	}, true);
	
}