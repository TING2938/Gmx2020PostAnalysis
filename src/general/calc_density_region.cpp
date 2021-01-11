/*  Author     : TING
 *  Date       : 2019/05/05
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate paired and unpair number along z axis in channel simulation.
 */

#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

};

gmx_main(temp)
{
	Handle hd(argc, argv);

	hd.ngrps = 1; /* grp1 and grp2, anion and cation. */

	/* add some user-defined pargs. */
	real dbin = 0.002;
	real lowPos = 0; // nm;
	real upPos = 30; // nm;
	int dim = 2;

	real lowPos2 = -2;
	real upPos2 = 20;
	int dim2 = 0;

	hd.pa = {
		{ "-dbin", FALSE, etREAL,  {&dbin}, "dbins."},
		{ "-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" },
		{ "-dim", FALSE, etINT, {&dim}, "dim"},
		{ "-up2", FALSE, etREAL, {&upPos2}, "up2 position of region of molecule/ion (nm)" },
		{ "-low2", FALSE, etREAL, {&lowPos2}, "low2 position of region of molecule/ion (nm)" },
		{ "-dim2", FALSE, etINT, {&dim2}, "dim2"}
	};

	hd.fnm = {
		{ efXVG, "-o", "density_region", ffWRITE }
	};

	hd.init();

	hd.readFirstFrame();

	auto posc1 = hd.initPosc(0);
	int nbin = (upPos - lowPos) / dbin;
	itp::vecd dens(nbin);
	dens.fill(0);
	int meZ;

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(posc1, 0);

		for (int i = 0; i != posc1.rows(); ++i)
		{
			if (lowPos <= posc1(i, dim) && posc1(i, dim) < upPos &&
				lowPos2 <= posc1(i, dim2) && posc1(i, dim2) < upPos2)
			{
				meZ = (posc1(i, dim) - lowPos) / dbin;
				dens[meZ] ++;
			}
		}

	} while (hd.readNextFrame());

	auto file = hd.openWrite(hd.get_opt2fn("-o"));
	fmt::print(file, "# number density (#)");
	for (int i = 0; i != nbin; ++i)
	{
		fmt::print(file, "{:10.4f} {:10d}\n", i * dbin + lowPos, dens[i]);
	}

	return 0;
}
