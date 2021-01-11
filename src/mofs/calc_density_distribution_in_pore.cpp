/*  Author     : TING
 *  Date       : 2019/06/05
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>

gmx_main(temp)
{
	itp::GmxHandle hd(argc, argv);

	// add some user-defined pargs.
	int nbin = 100;
	real lowPos = 7.761; // nm;
	real upPos = 12.239; // nm;
	real centerX = 2.829;
	real centerY = 2.056;
	real Rmof = 0.927;

	hd.pa = {
		{ "-nbin", FALSE, etINT,  {&nbin}, "nbins."},
		{ "-upPos", FALSE, etREAL, {&upPos}, "up position of region of pore (nm)" },
		{ "-lowPos", FALSE, etREAL, {&lowPos}, "low position of region of pore (nm)" },
		{ "-centerX", FALSE, etREAL, {&centerX}, "center X position of mofs pore (nm)" },
		{ "-centerY", FALSE, etREAL, {&centerY}, "center Y position of mofs pore (nm)" },
		{ "-Rmof", FALSE, etREAL, {&Rmof}, "Radius of mofs pore (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "numberDensity", ffWRITE }
	};
	hd.init();
	hd.readFirstFrame();

	auto posc = hd.initPosc(0);
	itp::vecd density(nbin);
	density.fill(0);

	double Rsquare = double(Rmof) * Rmof;
	double rsquare;
	double dbin = (upPos - lowPos) / nbin;
	size_t meZ = 0;

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(posc, 0);

		for (size_t i = 0; i != posc.rows(); ++i)
		{
			if (lowPos <= posc(i, ZZ) && posc(i, ZZ) < upPos)
			{
				rsquare = std::pow(posc(i, XX) - centerX, 2) + std::pow(posc(i, YY) - centerY, 2);
				if (rsquare <= Rsquare)
				{
					meZ = (size_t)((posc(i, ZZ) - lowPos) / dbin);
					density[meZ]++;
				}
			}
		}

	} while (hd.readNextFrame());

	double dvol = dbin * itp::pi * Rmof * Rmof;
	density /= dvol;

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	for (size_t i = 0; i != nbin; ++i)
	{
		fprintf(ofile, "%8.4e %8.4e\n", lowPos + dbin / 2 + i * dbin, density[i] / hd.nframe);
	}

	return 0;
}
