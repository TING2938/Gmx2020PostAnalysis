/*  Author     : TING
 *  Date       : 2019/06/05
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

public:
	itp::matd posc;
	itp::matd forcec;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

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
		{ efXVG, "-o", "pullforce", ffWRITE }
	};

	hd.flags = TRX_NEED_X | TRX_NEED_F;

	hd.init();
	hd.readFirstFrame();

	hd.posc.resize(hd.nmol[0], 3);
	hd.forcec.resize(hd.nmol[0], 3);
	
	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc, 0);
		hd.loadForceCenter(hd.forcec, 0);

		fmt::print(ofile, "{:12.8f} ", hd.fr->time);
		for (int i = 0; i < 3; i++)
		{
			fmt::print(ofile, "{:12.8f} ", hd.posc(0, i));
		}
		for (int i = 0; i < 3; i++)
		{
			fmt::print(ofile, "{:12.8f} ", hd.forcec(0, i));
		}
		fmt::print(ofile, "\n");

	} while (hd.readNextFrame());

	
	return 0;
}
