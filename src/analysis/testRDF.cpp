#include <itp/core>
#include <fstream>
#include <Eigen/Dense>

double safe_acos(double value)
{
	if (value <= -1.0)
	{
		return itp::pi;
	}
	else if (value >= 1.0)
	{
		return 0;
	}
	else
	{
		return std::acos(value);
	}
}

double calcVolume(double RMOF, double R, double L)
{
	double ret = 0;
	if (R <= RMOF - L)
	{
		ret = 4 * itp::pi * R * R * R / 3;
	}
	else if (R > RMOF + L)
	{
		double hmax = std::sqrt(R * R - (RMOF - L) * (RMOF - L));
		double hmin = std::sqrt(R * R - (RMOF + L) * (RMOF + L));
		double dh = 1e-5;
		double r, alpha, belta, cosAlpha, cosBelta;
		for (double h = hmin; h <= hmax; h += dh)
		{
			r = std::sqrt(R * R - h * h);
			cosAlpha = (r * r + L * L - RMOF * RMOF) / (2.0 * r * L);
			cosBelta = (L * L + RMOF * RMOF - r * r) / (2.0 * L * RMOF);
			alpha = safe_acos(cosAlpha);
			belta = safe_acos(cosBelta);
			ret += dh * (alpha * r * r + belta * RMOF * RMOF - L * RMOF * std::sin(belta));
		}
		ret += itp::pi / 3 * (3 * R - (R - hmax)) * (R - hmax) * (R - hmax);
		ret += hmin * itp::pi * RMOF * RMOF;
		ret *= 2;
	}
	else
	{
		double hmax = std::sqrt(R * R - (RMOF - L) * (RMOF - L));
		double dh = 1e-5;
		double r, alpha, belta, cosAlpha, cosBelta;
		for (double h = 0; h <= hmax; h += dh)
		{
			r = std::sqrt(R * R - h * h);
			cosAlpha = (r * r + L * L - RMOF * RMOF) / (2.0 * r * L);
			cosBelta = (L * L + RMOF * RMOF - r * r) / (2.0 * L * RMOF);
			alpha = safe_acos(cosAlpha);

			belta = safe_acos(cosBelta);

			ret += dh * (alpha * r * r + belta * RMOF * RMOF - L * RMOF * std::sin(belta));
		}
		ret += itp::pi / 3 * (3 * R - (R - hmax)) * (R - hmax) * (R - hmax);
		ret *= 2;
	}
	return ret;
}


int main()
{
	double upPos = 32.927;
	double lowPos = 27.073;

	double RMOF = 0.655;
	double L = 0.4;
	double dL = 0.001;

	int Lbin = RMOF / dL;
	double LMOF = upPos - lowPos;
	double Rmax = LMOF / 2;
	double dR = 0.002; // nm
	int Rbin = Rmax / dR;

	Eigen::ArrayXXd totVolume(Lbin, Rbin);
	totVolume.fill(0);
	fmt::print("{}\n", calcVolume(RMOF, 0.656, 0.001));
	

}