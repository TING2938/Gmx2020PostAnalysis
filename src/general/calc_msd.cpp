/*  Author     : TING
 *  Date       : 2020/01/09
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate msd.
 */

#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	void prepData(itp::matd& xcur, itp::matd& xprev, int grp)
	{
		double hbox[3];
		for (int m = 0; m != DIM; ++m)
		{
			hbox[m] = 0.5 * Lbox[m];
		}
		for (int i = 0; i != nmol[grp]; ++i)
		{
			for (int m = 0; m != DIM; ++m)
			{
				if (hbox[m] == 0)
				{
					continue;
				}
				while (xcur(i, m) - xprev(i, m) <= -hbox[m])
				{
					xcur(i, m) += Lbox[m];
				}
				while (xcur(i, m) - xprev(i, m) > hbox[m])
				{
					xcur(i, m) -= Lbox[m];
				}
			}
		}
	}

	itp::vecd calcMSD(const std::vector<itp::matd>& pos)
	{
		int nframe = pos.size();
		int halframe = nframe / 2;
		int nmol = pos[0].rows();

		itp::vecd msd(halframe);
		itp::veci count(halframe);

		msd.fill(0);
		count.fill(0);

		double dx, dy, dz;

		for (int i = 0; i != halframe; ++i)
		{
			for (int j = 0; j != halframe - i; ++j)
			{
				for (int k = 0; k != nmol; ++k)
				{
					dx = pos[i + j](k, XX) - pos[j](k, XX);
					dy = pos[i + j](k, YY) - pos[j](k, YY);
					dz = pos[i + j](k, ZZ) - pos[j](k, ZZ);
					msd[i] += dx * dx + dy + dy + dz + dz;
					count[i]++;
				}
			}
		}

		for (int i = 0; i != halframe; ++i)
		{
			if (count[i] != 0)
			{
				msd[i] /= count[i];
			}
		}
		return msd;
	}
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	/* add some user-defined pargs. */
	int nbin = 100;
	real lowPos = 0; // nm;
	real upPos = 30; // nm;
	const char* dimstr = "x";

	hd.pa = {
		{ "-d", FALSE, etSTR, {&dimstr}, "dim to calculate."},
		{ "-nbin", FALSE, etINT,  {&nbin}, "nbins."},
		{ "-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_msd", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();
	itp::matd xcur = hd.initPosc(0);   // the current coordinates 
	itp::matd xprev = hd.initPosc(0);  // the previous coordinates
	hd.loadPositionCenter(xprev, 0);

	std::vector<itp::matd> pos;

	double tmpTime = 0, dt;
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(xcur, 0);
		hd.prepData(xcur, xprev, 0);
		pos.push_back(xcur);

		xprev = xcur;
		dt = hd.fr->time - tmpTime;
		tmpTime = hd.fr->time;
	} while (hd.readNextFrame());

	auto msd = hd.calcMSD(pos);

	auto file = hd.openWrite(hd.get_ftp2fn(efXVG));

	for (int i = 0; i != msd.size(); ++i)
	{
		fmt::print(file, "{:8.4f} {:8.4f} \n", i * dt, msd[i]);
	}

	return 0;
}
