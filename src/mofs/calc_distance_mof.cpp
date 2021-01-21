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

	void getPosc(itp::vecd& posc, int grp);
public:
	itp::vecd posc1;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	/* add some user-defined pargs. */
	int nbin = 100;

	hd.pa = {
		{ "-nbin", FALSE, etINT,  {&nbin}, "nbins, useless"}
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_distance", ffWRITE }
	};
	
	hd.init();
	hd.readFirstFrame();
	
	hd.posc1.resize(3);

	auto file = hd.openWrite(hd.get_opt2fn("-o"));
	/* ---------------------------------------------------------------------------- */
	do {
		hd.getPosc(hd.posc1, 0);

		fmt::print(file, "{:8.4f} {:8.4f} {:8.4f} {:8.4f}\n", hd.fr->time, hd.posc1(XX), hd.posc1(YY), hd.posc1(ZZ));

	} while (hd.readNextFrame());

	return 0;
}

void Handle::getPosc(itp::vecd& posc, int grp)
{
	if (grp >= ngrps)
	{
		gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
	}

	int i, j, k, molN;

	double tmpPos, tmpCenter, halfLbox[3];

	itp::vecd atomCenter(ngx[grp]);
	for (i = 0; i != ngx[grp]; ++i)
	{
		atomCenter[i] = top->atoms.atom[index[grp][i]].m;
	}

	for (k = 0; k < 3; k++)
		halfLbox[k] = Lbox[k] / 2;

	molN = index[grp][0];
	for (k = 0; k != 3; k++)
	{
		tmpCenter = 0.0;  // initiate for each molecule
		for (j = 0; j != ngx[grp]; j++)
		{
			tmpPos = fr->x[molN + j][k];
			// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
			while (tmpPos - fr->x[molN][k] > halfLbox[k])
				tmpPos -= Lbox[k];
			while (tmpPos - fr->x[molN][k] < -halfLbox[k])
				tmpPos += Lbox[k];
			tmpCenter += tmpPos * atomCenter[j];
		}
		posc[k] = tmpCenter / atomCenter.sum();
	}
}
