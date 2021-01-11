/*  Author     : TING
 *  Date       : 2019/05/05
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate paired and unpair number along z axis in channel simulation.
 */
#include <itp/gmx>

int my_main(int argc, char* argv[])
{
	itp::GmxHandle hd(argc, argv);
	/* user-defined pargs. */ 
	int lateral = 2;
	int dim1 = 0;
	real lowPos1 = -2; // nm;
	real upPos1 = 30; // nm;
	real dim2 = 2;
	real lowPos2 = -2;
	real upPos2 = 30;

	hd.pa = {
		{ "-lateral", FALSE, etINT, {&lateral}, "lateral to calc. 0(x), 1(y), 2(z)"},
		{ "-dim1", FALSE, etINT, {&dim1}, "direction to calc. 0(x), 1(y), 2(z)"},
		{ "-up1", FALSE, etREAL, {&upPos1}, "up position of region of molecule/ion (nm)" },
		{ "-low1", FALSE, etREAL, {&lowPos1}, "low position of region of molecule/ion (nm)" },
		{ "-dim2", FALSE, etINT, {&dim2}, "direction to calc. 0(x), 1(y), 2(z)"},
		{ "-up2", FALSE, etREAL, {&upPos2}, "up position of region of molecule/ion (nm)" },
		{ "-low2", FALSE, etREAL, {&lowPos2}, "low position of region of molecule/ion (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "msd_lateral", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();

	matd posc_back = hd.initPosc(0);
	matd posc = hd.initPosc(0);
	int number;
	double res, r;

	hd.loadPositionCenter(posc_back, 0);
	hd.readNextFrame();
	auto file = hd.openWrite(hd.get_ftp2fn(efXVG));
	/* ---------------------------------------------------------------------------- */
	do {
		hd.loadPositionCenter(posc, 0);

		number = 0;
		res = 0;
		for (int i = 0; i != posc.nrow(); ++i) {
			if (lowPos1 < posc_back(i, dim1) && posc_back(i, dim1) < upPos1 && lowPos2 < posc_back(i, dim2) && posc_back(i, dim2) < upPos2) {
				r = 0;
				for (int m = 0; m != 3; ++m) {
					if (m != lateral) {
						r += std::pow(itp::periodicity(posc(i, m) - posc_back(i, m), hd.Lbox[m]), 2);
					}
				}
				res += r;
				number++;
			}
		}
		
		posc_back = posc;

		fprintf(file, "%12.4f %12.4f %10d\n", hd.fr->time, res / number, number);
		
	} while (hd.readNextFrame());

	/* ---------------------------------------------------------------------------- */
	
	/*********************************************************************************/
	return 0;
}



int main(int argc, char* argv[])
{
	return gmx_run_cmain(argc, argv, &my_main);
}
