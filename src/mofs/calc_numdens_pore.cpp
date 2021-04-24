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

	bool bInRegion(int num)
	{
		real tmp;
		for (int m = 0; m < 3; m++)
		{
			tmp = posc(num, m);
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

	hd.pa = {
		{ "-low", FALSE, etRVEC, {lowPos}, "low position" },
		{ "-up", FALSE, etRVEC, {upPos}, "up position" },
		{ "-com", FALSE, etINT, {&com}, "center of molecule, 0:mass, 1:geometry, 2:charge"}
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_numdens_pore", ffWRITE }
	};

	int grp = 0;
	hd.init();
	std::copy_n(lowPos, 3, hd.low);
	std::copy_n(upPos, 3, hd.up);

	hd.readFirstFrame();

	hd.posc = hd.initPosc(grp);

	double numdens = 0;
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc, grp, com);
		
		for (int i = 0; i < hd.nmol[grp]; i++)
		{
			if (hd.bInRegion(i))
			{
				numdens++;
			}
		}

	} while (hd.readNextFrame());

	numdens /= hd.nframe;

	auto file = hd.openWrite(hd.get_ftp2fn(efXVG));

	fmt::print(file, "{:8.4f}\n", numdens);


	return 0;
}
