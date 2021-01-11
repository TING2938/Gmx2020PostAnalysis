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

	void getPosc(Eigen::Array3d& posc, int grp);
public:
	Eigen::Array3d posc1, posc2;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	hd.ngrps = 2; /* number of group(s), grp1 and grp2, anion and cation. */

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
		{ efXVG, "-o", "calc_distance", ffWRITE }
	};
	
	hd.init(false);
	hd.readFirstFrame();
	
	hd.posc1.resize(3);
	hd.posc2.resize(3);

	auto file = hd.openWrite(hd.get_opt2fn("-o"));
	double d, dx, dy, dz;
	/* ---------------------------------------------------------------------------- */
	do {
		hd.getPosc(hd.posc1, 0);
		hd.getPosc(hd.posc2, 1);
		
		dx = hd.periodicity(hd.posc1(XX) - hd.posc2(XX), hd.Lbox[XX]);
		dy = hd.periodicity(hd.posc1(YY) - hd.posc2(YY), hd.Lbox[YY]);
		dz = hd.periodicity(hd.posc1(ZZ) - hd.posc2(ZZ), hd.Lbox[ZZ]);

		d = std::sqrt(dx * dx + dy * dy + dz * dz);
		fmt::print(file, "{:8.4f} {:8.4f}\n", hd.fr->time, d);


	} while (hd.readNextFrame());
	
	fclose(file);

	return 0;
}

void Handle::getPosc(Eigen::Array3d& posc, int grp)
{
	if (grp >= ngrps)
	{
		gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
	}

	int i, j, k, molN;

	double tmpPos, tmpCenter, halfLbox[3];

	Eigen::ArrayXd atomCenter(ngx[grp]);
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
