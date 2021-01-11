/*  Author     : TING
 *  Date       : 2019/06/04
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
		{ efXVG, "-o", "rdf_in_pore", ffWRITE }
	};
	hd.init(); 
	hd.readFirstFrame(); 

	auto posc = hd.initPosc(0);
	itp::vecd rdf(nbin);
	rdf.fill(0);

	double Rsquare = Rmof * Rmof;
	double rsquare;
	double dbin = Rmof / nbin;
	size_t meZ = 0;

	/* ---------------------------------------------------------------------------- */
	do { 
		hd.loadPositionCenter(posc, 0);

		for (size_t i = 0; i != posc.rows(); ++i) {
			if (lowPos <= posc(i, ZZ) && posc(i, ZZ) < upPos) {
				rsquare = std::pow(posc(i, XX) - centerX, 2) + std::pow(posc(i, YY) - centerY, 2);
				if (rsquare <= Rsquare) {
					meZ = (size_t)(std::sqrt(rsquare) / dbin);
					rdf[meZ]++;
				}
			}
		}  

	} while (hd.readNextFrame()); 

	itp::vecd volumn(nbin);
	volumn.fill(0);
	for (size_t i = 0; i != nbin; ++i) {
		volumn[i] = dbin * (i + 1) * dbin * (i + 1);
	}

	for (size_t i = nbin - 1; i != 0; --i) {
		volumn[i] = volumn[i] - volumn[i - 1];
	}

	volumn *= (itp::pi * (upPos - lowPos));
	rdf /= volumn;

	auto ofile = hd.openWrite(hd.get_opt2fn("-o")); 
	for (size_t i = 0; i != nbin; ++i) {
		fprintf(ofile, "%8.4e %8.4e\n", dbin / 2 + i * dbin, rdf[i] / hd.nframe);
	}

	return 0;
}
