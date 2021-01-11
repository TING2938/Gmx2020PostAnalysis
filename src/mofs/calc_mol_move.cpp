/*  Author     : TING
 *  Date       : 2019/06/05
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>

vector<int> autoACF(const vector<vector<int>>& data);

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
		{ efXVG, "-o", "mol_move", ffWRITE }
	};
	hd.init();
	hd.readFirstFrame();

	auto posc = hd.initPosc(0);

	vector<int> region(hd.nmol[0]);
	vector<int> regionPre(hd.nmol[0]);
	vector<double> cnt(6);
	std::fill(cnt.begin(), cnt.end(), 0.0);
	double Rsquare;

	gmx_bool bFirstFrame = 1;

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(posc, 0);

		regionPre = region;
		std::fill(region.begin(), region.end(), 0);

		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (lowPos <= posc(i, dim) && posc(i, dim) < upPos)
			{
				Rsquare = std::pow(posc(i, XX) - centerX, 2) + std::pow(posc(i, YY) - centerY, 2);
				if (Rsquare < Cr * Cr)
				{
					region[i] = 1;
				}
				else if (Rsquare <= CR * CR)
				{
					region[i] = 2;
				}
			}
		}

		if (bFirstFrame)
		{
			bFirstFrame = 0;
		}
		else
		{
			for (int i = 0; i < hd.nmol[0]; i++)
			{
				if (regionPre[i] == 0 && region[i] == 1)
					cnt[0]++;
				if (regionPre[i] == 1 && region[i] == 0)
					cnt[1]++;
				if (regionPre[i] == 0 && region[i] == 2)
					cnt[2]++;
				if (regionPre[i] == 2 && region[i] == 0)
					cnt[3]++;
				if (regionPre[i] == 1 && region[i] == 2)
					cnt[4]++;
				if (regionPre[i] == 2 && region[i] == 1)
					cnt[5]++;
			}
		}
	} while (hd.readNextFrame());

	for (auto&& c : cnt)
	{
		c /= (hd.nframe - 1);
	}

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	fmt::print(ofile, "# o-c  c-o  o-s  s-o  c-s  s-c\n");
	for (auto&& c : cnt)
	{
		fmt::print(ofile, "{:8.4f}\t", c);
	}
	fmt::print(ofile, "\n");

	return 0;
}

inline bool contain(const vector<int>& vec, int val)
{
	return std::find(vec.begin(), vec.end(), val) != vec.end();
}

vector<int> autoACF(const vector<vector<int>>& data)
{
	std::vector<int> acf;
	int nframe = data.size();
	int halframe = nframe / 2;
	acf.resize(halframe);
	std::fill(acf.begin(), acf.end(), 0);

	for (int i = 0; i < halframe; i++)
	{
		for (int m = 0; m < data[i].size(); m++)
		{
			for (int j = 0; j < halframe; j++)
			{
				if (contain(data[i + j], data[i][m]))
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
	return acf;
}
