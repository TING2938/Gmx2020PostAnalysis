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
	int grp = 0;

	/* add some user-defined pargs. */
	int com = 0;  // center of molecule, 0:mass, 1:geometry, 2:charge
	int nbin = 100;
	int dim = 2; // dim, 0(x), 1(y), 2(z);
	real lowPos = 0; // nm;
	real upPos = 30; // nm;
	
	hd.pa = {
		{ "-com", FALSE, etINT, {&com}, "center of molecule, 0:mass, 1:geometry, 2:charge"},
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

	hd.pos1.resize(hd.nmol[grp], hd.napm[grp]);
	hd.posc1.resize(hd.nmol[grp], 3);

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPosition(hd.pos1, grp);
		hd.loadPositionCenter(hd.posc1, grp, com);

	} while (hd.readNextFrame());

	return 0;
}
