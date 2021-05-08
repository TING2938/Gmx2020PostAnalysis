/*  Author     : TING
 *  Date       : 2019/09/09
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculation ecacf
 */
#include <itp/gmx>

class Handle : public itp::GmxHandleFull
{
public:
	using GmxHandleFull::GmxHandleFull;

	void rmPBC(int grp)
	{
		#pragma omp parallel for num_threads(nthreads)
		for (int j = 0; j < nmol[grp]; j++)
		{
			for (int m = 0; m < 3; m++)
			{
				while (posc(j, m) - posc_pre(j, m) > Lbox[m] / 2)
					posc(j, m) -= Lbox[m];
				while (posc(j, m) - posc_pre(j, m) < -Lbox[m] / 2)
					posc(j, m) += Lbox[m];
			}
		}
	}

	bool bInRegion(int nm)
	{
		double tmp;
		for (int m = 0; m < 3; m++)
		{
			tmp = posc(nm, m);
			while (tmp > Lbox[m])
				tmp -= Lbox[m];
			while (tmp < 0)
				tmp += Lbox[m];
			if (!(low[m] <= tmp && tmp < up[m]))
			{
				return false;
			}
		}
		return true;
	}

public:

	Eigen::ArrayXd cacf;
	std::vector<std::array<double, 3>> cacf1;

	itp::matd posc;
	itp::matd posc_pre;
	itp::matd velc;

	rvec low;
	rvec up;
};

gmx_main(temp)
{
	double kb = 1.38064852e-23; // J/K

	Handle hd(argc, argv);

	hd.flags = TRX_READ_X | TRX_READ_V;
	hd.ngrps = 1;
	rvec lowPos = { 0, 0, 0 };
	rvec upPos = { 30, 30, 30 };
	int com = 0;  // center of molecule, 0:mass, 1:geometry, 2:charge
	int Type = 6; // 0:X, 1:Y, 2:Z, 3:XY, 4:YZ, 5:XZ, 6:XYZ 
	real temperature = 300; // K
	real vol = 10; // nm^3
	std::vector<std::vector<int>> CalcDim = {
		{XX}, {YY}, {ZZ},
		{XX, YY}, {YY, ZZ}, {XX, ZZ},
		{XX, YY, ZZ},
	};

	hd.pa = {
		{"-temp", FALSE, etREAL, {&temperature}, "system temperature (K)"},
		{"-vol", FALSE, etREAL, {&vol}, "volume (nm^3)"},
		{ "-low", FALSE, etRVEC, { lowPos }, "low position" },
		{ "-up", FALSE, etRVEC, {upPos}, "up position" },
		{ "-com", FALSE, etINT, {&com}, "center of molecule, 0:mass, 1:geometry, 2:charge" },
		{ "-type", FALSE, etINT, {&Type}, "type to calculate. 0:X, 1:Y, 2:Z, 3:XY, 4:YZ, 5:XZ, 6:XYZ" }
	};

	hd.fnm = {
		{ efXVG, "-o", "IR", ffWRITE }
	};

	hd.init();
	std::copy_n(lowPos, 3, hd.low);
	std::copy_n(upPos, 3, hd.up);
	hd.readFirstFrame();

	hd.initPosc(hd.posc, 0);
	hd.initPosc(hd.posc_pre, 0);
	hd.initPosc(hd.velc, 0);
	hd.loadPositionCenter(hd.posc_pre, 0, com);

	double q = 0;
	double tmp;
	std::array<double, 3> qv;
	
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc, 0, com);
		hd.loadVelocityCenter(hd.velc, 0, com);

		hd.rmPBC(0);
		hd.posc_pre = hd.posc;

		std::fill_n(qv.begin(), 3, 0);

		for (int j = 0; j < hd.nmol[0]; j++) // for each molecule
		{
			if (hd.bInRegion(j))
			{
				q = hd.totCharge[0][j];
				for (int m = 0; m < 3; m++)
				{
					qv[m] += q * hd.velc(j, m);
				}
			}
		}
		hd.cacf1.push_back(qv);
	
	} while (hd.readNextFrame());

	real halframe = hd.nframe / 2;
	hd.cacf.resize(halframe);
	hd.cacf.fill(0);

	for (int i = 0; i < halframe; i++)
	{
		#pragma omp parallel for num_threads(nthreads)
		for (int j = 0; j < halframe; j++)
		{
			for (int m = 0; m < CalcDim[Type].size(); m++)
			{
				hd.cacf[j] += hd.cacf1[i][CalcDim[Type][m]] * hd.cacf1[i + j][CalcDim[Type][m]];
			}
		}
	}

	auto halfIndex = hd.nframe / 2;
	Eigen::ArrayXd time(halfIndex);
	for (int i = 0; i != halfIndex; i++)
	{
		hd.cacf[i] /= (halfIndex);
		time[i] = i * hd.dt;
	}
	
	Eigen::ArrayXd inteCacf = itp::cumtrapz(time, hd.cacf);
	Eigen::ArrayXd conductivity = inteCacf * (1 / (3 * vol * kb * temperature) * 2.56697e-17); // S/m

	auto file = hd.openWrite(hd.get_ftp2fn(efXVG));
	fmt::print(file, "# System volume is {:5.3f} nm^3\n", vol);
	fmt::print(file, "# System temperature is {:5.3f} K\n", temperature);
	fmt::print(file, "# time  cacf  conductivity (ps, e^2 nm^2 ps^-2, S/m)\n");

	for (int i = 0; i < halfIndex; i++)
	{
		fmt::print(file, "{:12.5f}  {:12.5f} {:12.5f}\n", time[i] , hd.cacf[i], conductivity[i]);
	}
	fclose(file);

	return 0;
}

