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
	itp::boxd pos1;
	itp::matd posc1;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	hd.ngrps = 1; /* number of group(s), grp1 and grp2, anion and cation. */

	/* add some user-defined pargs. */
	int nbin = 100;
	int dim = 2; // 0(x), 1(y), 2(z);
	real lowPos = 0; // nm;
	real upPos = 30; // nm;
	const char* str = "x";

	hd.pa = {
		{ "-str", FALSE, etSTR, {&str}, "str"},
		{ "-nbin", FALSE, etINT,  {&nbin}, "nbins."},
		{ "-dim", FALSE, etINT, {&dim}, "dim, 0(x), 1(y), 2(z)"},
		{ "-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_pair", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();

	hd.pos1.resize(hd.nmol[0], hd.napm[0]);
	hd.posc1.resize(hd.nmol[0], 3);

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPosition(hd.pos1, 0);
		hd.loadPositionCenter(hd.posc1, 0);

	} while (hd.readNextFrame());

	return 0;
}
