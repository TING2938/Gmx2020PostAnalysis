/*  Author     : TING
 *  Date       : 2019/09/09
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculation vacf
 */
#include <itp/gmx>
 
struct Handle : public itp::GmxHandle
{
	using GmxHandle::GmxHandle;
};

gmx_main(temp)
{

	double k = 0.999645; // kJ/mol;
	Handle hd(argc, argv);

	hd.flags = TRX_READ_X | TRX_READ_V;
	hd.ngrps = 1;

	int lateral = 2;
	int dim1 = 0;
	real lowPos1 = -2; // nm;
	real upPos1 = 30; // nm;
	int dim2 = 2;
	real lowPos2 = -2;
	real upPos2 = 30;

	hd.pa = {
		{ "-lateral", FALSE, etINT, {&lateral}, "lateral to calc. 0(x), 1(y), 2(z)"},
		{ "-d1", FALSE, etINT, {&dim1}, "direction to calc. 0(x), 1(y), 2(z)"},
		{ "-up1", FALSE, etREAL, {&upPos1}, "up position of region of molecule/ion (nm)" },
		{ "-low1", FALSE, etREAL, {&lowPos1}, "low position of region of molecule/ion (nm)" },
		{ "-d2", FALSE, etINT, {&dim2}, "direction to calc. 0(x), 1(y), 2(z)"},
		{ "-up2", FALSE, etREAL, {&upPos2}, "up position of region of molecule/ion (nm)" },
		{ "-low2", FALSE, etREAL, {&lowPos2}, "low position of region of molecule/ion (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "kinetic", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();
	auto mass = hd.mass[0][0].sum();

	auto velc = hd.initPosc(0);
	auto posc = hd.initPosc(0);
	
	double totEne = 0.0, ene = 0.0;
	int num = 0;
	auto file = hd.openWrite(hd.get_ftp2fn(efXVG));
	fmt::print(file, "# time(ps)  Energy(kJ/mol)\n");
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadVelocityCenter(velc, 0);
		hd.loadPositionCenter(posc, 0);

		ene = 0.0;
		num = 0;

		for (int i = 0; i != posc.nrow(); ++i)
		{
			if (lowPos1 < posc[i][dim1] && posc[i][dim1] < upPos1 && lowPos2 < posc[i][dim2] && posc[i][dim2] < upPos2)
			{
				for (int m = 0; m != 3; ++m)
				{
					if (m != lateral)
					{
						ene += std::pow(velc[i][m], 2);
					}
				}
				num++;
			}
		} 

		fmt::print(file, "{:.5f} {:8.3f}\n", hd.fr->time, 0.5 * k * mass * ene / num); 

	} while (hd.readNextFrame());

	fclose(file);

	return 0;
}

