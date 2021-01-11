/*  Author     : TING
 *  Date       : 2019/06/04
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

public:
	itp::matd posc;
};


gmx_main(temp)
{
	Handle hd(argc, argv);

	// add some user-defined pargs.
	real x1 = 0;
	real x2 = 0;
	real y1 = 0;
	real y2 = 0;
	real z1 = 0;
	real z2 = 0;
	real centerX = 2.829;
	real centerY = 2.056;
	int dim = 2;
	real dr = 0.01;
	real Cr = 0;
	real CR = 0.5;

	hd.pa = {
		{ "-x1", FALSE, etREAL, {&x1}, "low x position of region of pore (nm)" },
		{ "-x2", FALSE, etREAL, {&x2}, "up x position of region of pore (nm)" },
		{ "-y1", FALSE, etREAL, {&y1}, "low y position of region of pore (nm)" },
		{ "-y2", FALSE, etREAL, {&y2}, "up y position of region of pore (nm)" },
		{ "-z1", FALSE, etREAL, {&z1}, "low z position of region of pore (nm)" },
		{ "-z2", FALSE, etREAL, {&z2}, "up z position of region of pore (nm)" },
		{ "-centerX", FALSE, etREAL, {&centerX}, "center X position of mofs pore (nm)" },
		{ "-centerY", FALSE, etREAL, {&centerY}, "center Y position of mofs pore (nm)" },
		{"-dim", FALSE, etINT, {&dim}, "dim to analysis. 0, x; 1, y; 2, z"},
		{ "-dr", FALSE, etREAL, {&dr}, "dr of bin along radius (nm)" },
		{ "-Cr", FALSE, etREAL, {&Cr}, "min R region of pore (nm)" },
		{ "-CR", FALSE, etREAL, {&CR}, "max R region of pore (nm)" },
	};

	hd.fnm = {
		{ efXVG, "-o", "integration_cdf_in_pore", ffWRITE }
	};
	hd.init();
	hd.readFirstFrame();

	hd.posc.resize(hd.nmol[0], 3);
	int nbin = (CR - Cr) / dr;
	itp::vecd rdf(nbin);
	rdf.fill(0);

	itp::vecd volumn(nbin);
	volumn.fill(0);
	for (size_t i = 0; i != nbin; ++i)
	{
		volumn[i] = itp::pi * (std::pow(Cr + dr * (i + 1), 2) - std::pow(Cr + dr * i, 2)) * (z2 - z1);
	}

	auto file = hd.openWrite(hd.get_opt2fn("-o"));

	size_t meZ = 0;
	double r;
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc, 0);

		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (x1 <= hd.posc(i, XX) && hd.posc(i, XX) < x2 &&
				y1 <= hd.posc(i, YY) && hd.posc(i, YY) < y2 &&
				z1 <= hd.posc(i, ZZ) && hd.posc(i, ZZ) < z2)
			{
				r = std::sqrt(std::pow(hd.posc(i, XX) - centerX, 2) + std::pow(hd.posc(i, YY) - centerY, 2));
				meZ = (r - Cr) / dr;
				rdf[meZ]++;
			}
		}

	} while (hd.readNextFrame());

	rdf /= (hd.nframe * (z2-z1));
	
	auto cdf = itp::cumsum(rdf);

	for (int i = 0; i != nbin; ++i)
	{
		fmt::print(file, "{:8.4f} {:8.4f}\n", i*dr+dr/2, cdf[i]);
	}

	return 0;
}
