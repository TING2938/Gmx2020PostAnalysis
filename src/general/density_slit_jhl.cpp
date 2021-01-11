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

public:
	itp::matd posc1;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	hd.ngrps = 1; /* number of group(s), grp1 and grp2, anion and cation. */

	/* add some user-defined pargs. */
	int nbin1 = 100;
	int dim1 = 0; // 0(x), 1(y), 2(z);
	int dim2 = 2;
	real lowPos1 = 0; // nm;
	real upPos1 = 30; // nm;
	real upPos2 = 0;
	real lowPos2 = 10;
	real yy = 0;

	hd.pa = {
		{ "-nbin1", FALSE, etINT,  {&nbin1}, "nbins."},
		{ "-dim1", FALSE, etINT, {&dim1}, "dim, 0(x), 1(y), 2(z)"},
		{ "-dim2", FALSE, etINT, {&dim2}, "dim, 0(x), 1(y), 2(z)"},
		{ "-up1", FALSE, etREAL, {&upPos1}, "up position of region of molecule/ion (nm)" },
		{ "-low1", FALSE, etREAL, {&lowPos1}, "low position of region of molecule/ion (nm)" },
		{ "-up2", FALSE, etREAL, {&upPos2}, "up position of region of molecule/ion (nm)" },
		{ "-low2", FALSE, etREAL, {&lowPos2}, "low position of region of molecule/ion (nm)" },
		{ "-yy", FALSE, etREAL, {&yy}, "y length (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "density_slit", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();

	hd.posc1.resize(hd.nmol[0], 3);
	int me = 0;
	double dbin = (upPos1 - lowPos1) / nbin1;
	std::vector<double> density(nbin1, 0);
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc1, 0);
		for (int i = 0; i < hd.nmol[0]; i++)
		{
			if (lowPos1 <= hd.posc1(i, dim1) && hd.posc1(i, dim1) <= upPos1
				&& lowPos2 <= hd.posc1(i, dim2) && hd.posc1(i, dim2) <= upPos2)
			{
				me = (hd.posc1(i, dim1) - lowPos1) / dbin;
				density[me]++;
			}
		}

	} while (hd.readNextFrame());

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));

	double dVolume = (upPos2 - lowPos2) * dbin * yy;
	for (auto&& d : density)
	{
		d /= (hd.nframe * dVolume); // #/nm^3
	}

	for (int i = 0; i < nbin1; i++)
	{
		fmt::print(ofile, "{:15.8f} {:15.8f}\n", dbin/2 + dbin*i, density[i]);
	}
	return 0;
}
