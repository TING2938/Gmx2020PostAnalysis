/*  Author     : TING
 *  Date       : 2020/01/09
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate msd.
 */

#include <itp/gmx>

// 3D array
using array3d = std::array<double, 3>;

class Handle : public itp::GmxHandleFull
{
public:
	using GmxHandleFull::GmxHandleFull;

	void rmPBC(int grp)
	{
		for (int i = 1; i < nframe; i++)
		{
			for (int j = 0; j < nmol[grp]; j++)
			{
				for (int m = 0; m < 3; m++)
				{
					while (allPosc[i](j, m) - allPosc[i - 1](j, m) > Lbox[m] / 2)
						allPosc[i](j, m) -= Lbox[m];
					while (allPosc[i](j, m) - allPosc[i - 1](j, m) < -Lbox[m] / 2)
						allPosc[i](j, m) += Lbox[m];
				}
			}
		}
	}

	itp::vecd calcQMSD(const std::vector<int>& calcDim, int grp)
	{
		int nframe = allPosc.size();
		int halframe = nframe / 2;

		itp::matd qmsd(halframe, calcDim.size());
		itp::vecd count(halframe);
		qmsd.fill(0);
		count.fill(0);

		for (int i = 0; i < halframe; i++) // for each halframe
		{
			for (int j = 0; j < nmol[grp]; j++) // for each molecule
			{
				if (bInRegion(i, j))
				{
					for (int k = 0; k < halframe; k++)
					{
						if (bInRegion(i + k, j))
						{
							for (int m = 0; m < calcDim.size(); m++)
							{
								qmsd(k, m) += totCharge[grp][j] * (allPosc[i + k](j, calcDim[m]) - allPosc[i](j, calcDim[m]));
							}
							count[k]++;
						}
						else
						{
							break;
						}
					}
				}
			}
		}

		itp::vecd outputQMSD(halframe);
		outputQMSD.fill(0);

		for (int i = 0; i != halframe; ++i)
		{
			for (int m = 0; m < calcDim.size(); m++)
			{
				outputQMSD[i] += std::pow(qmsd(i, m), 2);
			}
			if (count[i] != 0)
				outputQMSD[i] /= count[i];
		}
		return outputQMSD;
	}

	bool bInRegion(int frame, int nm)
	{
		double tmp;
		for (int m = 0; m < 3; m++)
		{
			tmp = allPosc[frame](nm, m);
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
	itp::matd posc;
	std::vector<itp::matd> allPosc;
	rvec low;
	rvec up;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	/* add some user-defined pargs. */
	rvec lowPos = { 0, 0, 0 };
	rvec upPos = { 30, 30, 30 };
	int com = 0;  // center of molecule, 0:mass, 1:geometry, 2:charge
	int qmsdType = 6; // 0:X, 1:Y, 2:Z, 3:XY, 4:YZ, 5:XZ, 6:XYZ 
	std::vector<std::vector<int>> qmsdCalcDim = {
		{XX}, {YY}, {ZZ},
		{XX, YY}, {YY, ZZ}, {XX, ZZ},
		{XX, YY, ZZ},
	};

	hd.pa = {
		{ "-low", FALSE, etRVEC, {lowPos}, "low position" },
		{ "-up", FALSE, etRVEC, {upPos}, "up position" },
		{ "-com", FALSE, etINT, {&com}, "center of molecule, 0:mass, 1:geometry, 2:charge"},
		{ "-type", FALSE, etINT, {&qmsdType}, "type to calculate. 0:X, 1:Y, 2:Z, 3:XY, 4:YZ, 5:XZ, 6:XYZ"},
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_qmsd_Ft", ffWRITE }
	};

	int grp = 0;
	hd.init();
	std::copy_n(lowPos, 3, hd.low);
	std::copy_n(upPos, 3, hd.up);

	hd.readFirstFrame();

	hd.initPosc(hd.posc, grp);

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc, grp, com);
		hd.allPosc.push_back(hd.posc);
	} while (hd.readNextFrame());

	hd.rmPBC(grp);
	itp::vecd qmsd = hd.calcQMSD(qmsdCalcDim[qmsdType], grp);

	auto file = hd.openWrite(hd.get_ftp2fn(efXVG));

	for (int i = 0; i != qmsd.size(); ++i)
	{
		fmt::print(file, "{:8.4f} {:8.4f} \n", i * hd.dt, qmsd[i]);
	}

	return 0;
}
