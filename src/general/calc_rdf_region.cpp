/*  Author     : TING
 *  Date       : 2019/10/31
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate rdf of region
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

	hd.ngrps = 2; /* grp1 and grp2, anion and cation. */

	/* add some user-defined pargs. */
	int nbin = 5;
	real lowPos = -100; // nm;
	real upPos = 100; // nm;
	int dim = 0; // 0(x), 1(y), 2(z)

	int dim2 = 2;
	real lowPos2 = -100, upPos2 = 100;
	real Rm = 0;

	hd.pa = {
		{ "-nbin", FALSE, etINT,  {&nbin}, "nbins."},
		{ "-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" },
		{ "-dim", FALSE, etINT, {&dim}, "dim of up and low position, 0(x), 1(y), 2(z)"},

		{ "-dim2", FALSE, etINT, {&dim2}, "dim2"},
		{ "-low2", FALSE, etREAL, {&lowPos2}, "z1"},
		{ "-up2", FALSE, etREAL, {&upPos2}, "z2"},
		{ "-Rm", FALSE, etREAL, {&Rm}, "Rm"}
	};

	hd.fnm = {
		{ efXVG, "-o", "rdf_region", ffWRITE }
	};

	hd.init();

	hd.readFirstFrame();

	if (Rm < 1e-6)
	{
		Rm = std::min(hd.Lbox[0], hd.Lbox[1]);
		Rm = std::min(double(Rm), hd.Lbox[2]);
		Rm /= 2;
	}
	int Rbin = Rm / 0.002;
	double R, r;

	auto posc1 = hd.initPosc(0);
	auto posc2 = hd.initPosc(1);
	itp::matd rdf(nbin, Rbin);
	rdf.fill(0);
	int meR, meZ;

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(posc1, 0);
		hd.loadPositionCenter(posc2, 1);

		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (lowPos < posc1(i, dim) && posc1(i, dim) < upPos && lowPos2 < posc1(i, dim2) && posc1(i, dim2) < upPos2)
			{
				meZ = nbin * (posc1(i, dim) - lowPos) / (upPos - lowPos);

				for (int j = 0; j != hd.nmol[1]; ++j)
				{
					R = 0;
					for (int k = 0; k != 3; ++k)
					{
						r = hd.periodicity(posc1(i, k) - posc2(j, k), hd.Lbox[k]);
						R += r * r;
					}

					if (1e-6 < R && R < Rm * Rm)
					{
						meR = Rbin * std::sqrt(R) / Rm;
						rdf(meZ, meR)++;
					}
				}
			}
		}

	} while (hd.readNextFrame());


	itp::vecd volumn(Rbin);
	for (int i = 0; i != Rbin; ++i)
	{
		volumn[i] = std::pow(0.002 * (i + 1), 3);
	}

	for (int i = Rbin - 1; i != 0; --i)
	{
		volumn[i] = volumn[i] - volumn[i - 1];
	}

	rdf.rowwise() /= volumn.transpose();

	itp::vecd mean = rdf.rightCols(Rbin * 0.1).rowwise().mean();

	rdf.colwise() /= mean;

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));

	fmt::print(ofile, "# {} ", char('x' + dim));
	for (int i = 0; i != nbin + 1; ++i)
	{
		fmt::print(ofile, "{:8.4f} ", lowPos + i * (upPos - lowPos) / nbin);
	}
	fmt::print(ofile, "\n");

	for (int j = 0; j != Rbin; ++j)
	{
		fmt::print(ofile, "{:8.4f} ", j * 0.002);

		for (int i = 0; i != nbin; ++i)
		{
			fmt::print(ofile, "{:8.4f} ", rdf(i, j));
		}
		fmt::print(ofile, "\n");
	}
	return 0;
}
