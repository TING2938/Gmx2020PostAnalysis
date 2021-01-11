/*  Author     : TING
 *  Date       : 2019/06/05
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>

vector<double> getACF(const vector<vector<int>>& data, int type);

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
		{ efXVG, "-o", "inPoreNumAcf", ffWRITE }
	};
	hd.init();
	hd.readFirstFrame();

	auto posc = hd.initPosc(0);

	int nm = hd.nmol[0];
	vector<vector<int>> range;
	vector<int> tmp_range(nm);
	double Rsquare;
	

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(posc, 0);

		std::fill(tmp_range.begin(), tmp_range.end(), 0);
		
		for (int i = 0; i != nm; ++i)
		{
			if (lowPos <= posc(i, dim) && posc(i, dim) < upPos)
			{
				Rsquare = std::pow(posc(i, XX) - centerX, 2) + std::pow(posc(i, YY) - centerY, 2);
				if (Rsquare <= Cr * Cr)
				{
					tmp_range[i] = 1; // center region
				}
				else if (Rsquare <= CR * CR)
				{
					tmp_range[i] = 2; // surface region
				}
			}
		}
		range.push_back(tmp_range);

	} while (hd.readNextFrame());

	auto acf_center = getACF(range, 1);
	auto acf_surface = getACF(range, 2);

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	fmt::print(ofile, "# time   acf_center   acf_surface\n");
	for (int i = 0; i != acf_center.size(); i++)
	{
		fmt::print(ofile, "{:8.3f}\t{:8.5f}\t{:8.5f}\n", i * hd.dt, acf_center[i], acf_surface[i]);
	}
	return 0;
}

vector<double> getACF(const vector<vector<int>>& data, int type)
{

	int nframe = data.size();
	int nmol = data[0].size();
	vector<double> acf(nframe, 0.0);

	for (int i = 0; i < nframe; i++)
	{
		for (int m = 0; m < nmol; m++)
		{
			if (data[i][m] == type)
			{
				for (int j = 0; j < nframe - i; j++)
				{
					if (data[i + j][m] == type)
					{
						acf[j]++;
					}
					else
					{
						break;
					}
				}
			}
		}
	}
	for (int i = 0; i < nframe; i++)
	{
		acf[i] /= (nframe - i);
	}
	return acf;
}