/*  Author     : TING
 *  Date       : 2019/06/05
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>

gmx_main(temp)
{
	itp::GmxHandle hd(argc, argv);

	// add some user-defined pargs.
	real lowPos = 7.761; // nm;
	real upPos = 12.239; // nm;
	real dim = 2;
	real centerX = 2.829;
	real centerY = 2.056;
	real CR = 0.927;
	real Cr = 0.4;

	hd.pa = {
		{ "-upPos", FALSE, etREAL, {&upPos}, "up position of region of pore (nm)" },
		{ "-lowPos", FALSE, etREAL, {&lowPos}, "low position of region of pore (nm)" },
		{ "-dim", FALSE, etREAL, {&dim}, "dim, 0(x), 1(y), 2(z)" },
		{ "-centerX", FALSE, etREAL, {&centerX}, "center X position of mofs pore (nm)" },
		{ "-centerY", FALSE, etREAL, {&centerY}, "center Y position of mofs pore (nm)" },
		{ "-CR", FALSE, etREAL, {&CR}, "Radius of mofs pore (nm)" },
		{ "-Cr", FALSE, etREAL, {&Cr}, "Radius of pore central (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-center", "mol_index_center", ffWRITE },
		{ efXVG, "-surface", "mol_index_surface", ffWRITE }
	};
	hd.init();
	hd.readFirstFrame();

	auto posc = hd.initPosc(0);

	vector<int> region(hd.nmol[0]);
	vector<double> cnt(6);
	std::fill(cnt.begin(), cnt.end(), 0.0);
	double Rsquare;

	gmx_bool bFirstFrame = 1;

	auto fp_center = hd.openWrite(hd.get_opt2fn("-center"));
	auto fp_surface = hd.openWrite(hd.get_opt2fn("-surface"));

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(posc, 0);

		fmt::print(fp_center, "{:12.8f} ", hd.fr->time);
		fmt::print(fp_surface, "{:12.8f} ", hd.fr->time);

		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (lowPos <= posc(i, dim) && posc(i, dim) < upPos)
			{
				Rsquare = std::pow(posc(i, XX) - centerX, 2) + std::pow(posc(i, YY) - centerY, 2);
				if (Rsquare < Cr * Cr)
				{
					fmt::print(fp_center, "{:6d} ", i);
				}
				else if (Rsquare <= CR * CR)
				{
					fmt::print(fp_surface, "{:6d} ", i);
				}
			}
		}

		fmt::print(fp_center, "\n");
		fmt::print(fp_surface, "\n");
	} while (hd.readNextFrame());

	return 0;
}
