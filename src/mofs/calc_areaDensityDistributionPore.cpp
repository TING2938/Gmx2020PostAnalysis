/*  Author     : TING
 *  Date       : 2019/06/04
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>

gmx_main(temp)
{
	itp::GmxHandle hd(argc, argv);

	// add some user-defined pargs.
	int dim1 = 0;
	real lowPos1 = 1.66; // nm;
	real upPos1 = 3.996; // nm;

	int dim2 = 1;
	real lowPos2 = 1.133;
	real upPos2 = 3.073;

	int dim3 = 2;
	real lowPos3 = 7.761; // nm;
	real upPos3 = 12.239; // nm;

	real dbin = 0.002;

	hd.pa = {
		{ "-dim1", FALSE, etINT,  {&dim1}, "dim1."},
		{ "-up1", FALSE, etREAL, {&upPos1}, "up1 position of region of pore (nm)" },
		{ "-low1", FALSE, etREAL, {&lowPos1}, "low1 position of region of pore (nm)" },
		{ "-dim2", FALSE, etINT,  {&dim2}, "dim2."},
		{ "-up2", FALSE, etREAL, {&upPos2}, "up2 position of region of pore (nm)" },
		{ "-low2", FALSE, etREAL, {&lowPos2}, "low2 position of region of pore (nm)" },
		{ "-dim3", FALSE, etINT,  {&dim3}, "dim3. 7.761-12.239 or 27.76-32.239"},
		{ "-up3", FALSE, etREAL, {&upPos3}, "up3 position of region of pore (nm)" },
		{ "-low3", FALSE, etREAL, {&lowPos3}, "low3 position of region of pore (nm)" },
		{ "-dbin", FALSE, etREAL, {&dbin}, "dbin"}
	};

	hd.fnm = {
		{ efXVG, "-o", "areaDensity", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();

	auto posc = hd.initPosc(0);

	int nbin1 = (upPos1 - lowPos1) / dbin;
	int nbin2 = (upPos2 - lowPos2) / dbin;

	itp::matd dens(nbin1, nbin2);
	dens.fill(0);

	int me1, me2;

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(posc, 0);

		for (int i = 0; i != posc.rows(); ++i)
		{
			if (lowPos1 <= posc(i, dim1) && posc(i, dim1) < upPos1 &&
				lowPos2 <= posc(i, dim2) && posc(i, dim2) < upPos2 &&
				lowPos3 <= posc(i, dim3) && posc(i, dim3) < upPos3)
			{
				me1 = (posc(i, dim1) - lowPos1) / dbin;
				me2 = (posc(i, dim2) - lowPos2) / dbin;
				dens(me1, me2)++;
			}
		}
	} while (hd.readNextFrame());

	dens /= (dbin * dbin * (upPos3 - lowPos3) * hd.nframe);

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	fmt::print(ofile, "# #/nm^3\n");

	for (int i = 0; i != nbin1; ++i)
	{
		for (int j = 0; j != nbin2; ++j)
		{
			fmt::print(ofile, "{:8.4f} ", dens(i, j));
		}
		fmt::print(ofile, "\n");
	}

	return 0;
}
