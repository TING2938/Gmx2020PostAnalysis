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
	real centerX = 2.829;
	real centerY = 2.056;
	real Rmof = 0.927;

	hd.pa = {
		{ "-upPos", FALSE, etREAL, {&upPos}, "up position of region of pore (nm)" },
		{ "-lowPos", FALSE, etREAL, {&lowPos}, "low position of region of pore (nm)" },
		{ "-centerX", FALSE, etREAL, {&centerX}, "center X position of mofs pore (nm)" },
		{ "-centerY", FALSE, etREAL, {&centerY}, "center Y position of mofs pore (nm)" },
		{ "-Rmof", FALSE, etREAL, {&Rmof}, "Radius of mofs pore (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "DensityTime", ffWRITE }
	};
	hd.init();
	hd.readFirstFrame();

	auto posc = hd.initPosc(0);

	double Rsquare = double(Rmof) * Rmof;
	double rsquare;
	int density;
	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(posc, 0);

		density = 0;
		for (size_t i = 0; i != posc.rows(); ++i)
		{
			if (lowPos <= posc(i, ZZ) && posc(i, ZZ) < upPos)
			{
				rsquare = std::pow(posc(i, XX) - centerX, 2) + std::pow(posc(i, YY) - centerY, 2);
				if (rsquare <= Rsquare)
				{
					density++;
				}
			}
		}

		fmt::print(ofile, "{:8.4f} {:10d}\n", hd.fr->time, density);

	} while (hd.readNextFrame());

	fclose(ofile);

	return 0;
}
