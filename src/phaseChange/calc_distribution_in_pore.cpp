/*  Author     : TING
 *  Date       : 2019/06/04
 *  Email      : yeting2938@hust.edu.cn 
 */

#include <itp/gmx>

void getElectrodePosition(itp::GmxHandle& hd, int grp1, int grp2, double& lowPos, double& upPos,
	double& centerX, double& centerY, double& aX, double& aY, double factor);

gmx_main(temp)
{
	itp::GmxHandle hd(argc, argv);

	hd.ngrps = 3;

	// add some user-defined pargs.
	int nbin = 100;
	double centerX, centerY;
	double lowPos, upPos;
	double aX, aY;

	hd.pa = {
		{ "-nbin", FALSE, etINT,  {&nbin}, "nbins."}
	};


	hd.fnm = {
		{ efXVG, "-o", "numberDensity", ffWRITE }
	};
	hd.init();
	hd.readFirstFrame();

	getElectrodePosition(hd, 1, 2, lowPos, upPos, centerX, centerY, aX, aY, 0.8);

	int nm = hd.nmols[0]; // nr. of molecule
	matd posc(nm, 3);
	vecd density(nbin, 0);

	double dbin = (upPos - lowPos) / nbin;
	size_t meZ;

	/* ---------------------------------------------------------------------------- */
	do {
		hd.loadPositionCenter(posc, 0);

		for (size_t i = 0; i != nm; ++i) {
			if (lowPos <= posc(i, ZZ) && posc(i, ZZ) < upPos && std::abs(posc(i, XX) - centerX) <= aX / 2 && std::abs(posc(i, YY) - centerY) <= aY / 2) {
				meZ = (size_t)((posc(i, ZZ) - lowPos) / dbin);
				++density[meZ];
			}
		}
	} while (hd.readNextFrame());

	double dvol = dbin * aX * aY;
	density /= dvol;

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	for (size_t i = 0; i != nbin; ++i) {
		fprintf(ofile, "%8.4e %8.4e\n", lowPos + dbin / 2 + i * dbin, density[i] / hd.nframe);
	}

	return 0;
}

void getElectrodePosition(itp::GmxHandle & hd, int grp1, int grp2, double & lowPos, double & upPos,
	double & centerX, double & centerY, double & aX, double& aY, double factor)
{
	int nm1 = hd.nmols[grp1];
	int nm2 = hd.nmols[grp2];
	matd pos1(nm1, 3), pos2(nm2, 3);

	hd.loadPositionCenter(pos1, grp1);
	hd.loadPositionCenter(pos2, grp2);

	lowPos = pos1(0, ZZ);
	upPos = pos2(0, ZZ);
	if (lowPos > upPos)
		std::swap(lowPos, upPos);

	double x1 = pos1.min(XX);
	double x2 = pos1.max(XX);
	double y1 = pos1.min(YY);
	double y2 = pos1.max(YY);

	centerX = (x1 + x2) / 2;
	centerY = (y1 + y2) / 2;

	aX = (x2 - x1) * factor;
	aY = (y2 - y1) * factor;
}