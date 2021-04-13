/*  Author     : TING
 *  Date       : 2021/04/13
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate msd.
 */

#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	void rmPBC(std::vector<itp::matd>& allPos, int grp)
	{
		for (int i = 1; i < nframe; i++)
		{
			for (int j = 0; j < nmol[grp]; j++)
			{
				for (int m = 0; m < 3; m++)
				{
					while (allPos[i](j, m) - allPos[i - 1](j, m) > Lbox[m] / 2)
						allPos[i](j, m) -= Lbox[m];
					while (allPos[i](j, m) - allPos[i - 1](j, m) < -Lbox[m] / 2)
						allPos[i](j, m) += Lbox[m];
				}
			}
		}
	}

	int getMe(real pos, real lowPos, double dbin, int dim)
	{
		while (pos > Lbox[dim])
			pos -= Lbox[dim];
		while (pos < 0)
			pos += Lbox[dim];
		return int((pos - lowPos) / dbin);
	}

	// return: msd
	itp::matd calcMSD(const std::vector<itp::matd>& allPos, const std::vector<int>& calcDim, 
		int grp, real lowPos, real upPos, int nbin, int dim)
	{
		int halframe = nframe / 2;
		double dbin = (upPos - lowPos) / nbin;
		itp::matd msd(halframe, nbin);
		itp::veci count(nbin);
		msd.fill(0);
		count.fill(0);
		int me;

		for (int k = 0; k < nmol[grp]; k++)
		{
			for (int i = 0; i < halframe; i++)
			{
				// determine which bin the position belong to
				me = getMe(allPos[i](k, dim), lowPos, dbin, dim);
				if (0 <= me && me < nbin)
				{
					count[me]++;
					for (int j = 0; j < halframe; j++)
					{
						for (auto&& m : calcDim)
						{
							msd(j, me) += std::pow(allPos[i + j](k, m) - allPos[i](k, m), 2);
						}
					}
				}
			}
		}

		for (int i = 0; i < halframe; i++)
		{
			for (int j = 0; j < nbin; j++)
			{
				if (count[j] > 0.0)
				{
					msd(i, j) /= count[j];
				}
			}
		}
		return msd;
	}

};

gmx_main(temp)
{
	Handle hd(argc, argv);

	/* add some user-defined pargs. */
	int dim = 2;
	int nbin = 100;
	real lowPos = 0; // nm;
	real upPos = 30; // nm;
	int msdType = 6; // 0:X, 1:Y, 2:Z, 3:XY, 4:YZ, 5:XZ, 6:XYZ 
	std::vector<std::vector<int>> msdCalcDim = {
		{0}, {1}, {2},
		{0, 1}, {1, 2}, {0, 2},
		{0, 1, 2},
	};

	hd.pa = {
		{ "-type", FALSE, etINT, {&msdType}, "type to calculate. 0:X, 1:Y, 2:Z, 3:XY, 4:YZ, 5:XZ, 6:XYZ"},
		{ "-dim", FALSE, etINT, {&dim}, "dim of selected region"},
		{ "-nbin", FALSE, etINT,  {&nbin}, "nbins."},
		{ "-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "msd_bin", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();

	int grp = 0;
	itp::matd posc = hd.initPosc(grp);
	std::vector<itp::matd> allPosc;

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(posc, 0);
		allPosc.push_back(posc);
	} while (hd.readNextFrame());

	hd.rmPBC(allPosc, grp);

	itp::matd msd = hd.calcMSD(allPosc, msdCalcDim[msdType], grp, lowPos, upPos, nbin, dim);

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	for (int i = 0; i < msd.rows(); i++)
	{
		fmt::print(ofile, "{:8.3f}\t", i * hd.dt);
		for (int j = 0; j < msd.cols(); j++)
		{
			fmt::print(ofile, "{:10.6f}\t", msd(i, j));
		}
		fmt::print(ofile, "\n");
	}

	return 0;
}
