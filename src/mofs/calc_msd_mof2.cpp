/*  Author     : TING
 *  Date       : 2020/01/09
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate msd.
 */

#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	/**
	 * @brief 计算每个原子在区域内能持续过久
	*/
	std::vector<std::vector<int>> calcDuration(const std::vector<std::vector<int>>& index)
	{
		int nframe = index.size();
		int halframe = nframe / 2;
		std::vector<std::vector<int>> dura(halframe);
		std::vector<int> tmpdura;

		for (int i = 0; i != halframe; ++i)
		{
			tmpdura.resize(index[i].size());
			std::fill(tmpdura.begin(), tmpdura.end(), 0);

			for (int j = 0; j != index[i].size(); ++j)
			{
				for (int k = 0; k != halframe; ++k)
				{
					if (itp::contains(index[i + k], index[i][j]))
					{
						tmpdura[j]++;
					}
					else
					{
						break;
					}
				}
			}
			dura[i].swap(tmpdura);
		}
		return dura;
	}

	itp::vecd calcMSD(const std::vector<std::vector<double>>& pos,
		const std::vector<std::vector<int>>& index, const std::vector<std::vector<int>>& duration)
	{
		int nframe = pos.size();
		int halframe = nframe / 2;

		itp::vecd msd = itp::vecd::Zero(halframe);
		itp::veci count = itp::veci::Zero(halframe);

		for (int i = 0; i != halframe; ++i)
		{
			for (int j = 0; j != index[i].size(); ++j)
			{
				for (int k = 0; k != duration[i][j]; ++k)
				{
					double pos1 = pos[i][j];
					double pos2 = pos[i + k][std::find(index[i + k].begin(), index[i + k].end(), index[i][j]) - index[i + k].begin()];
					msd[k] += ((pos1 - pos2) * (pos1 - pos2));
					count[k]++;
				}
			}
		}

		fmt::print("\n\nnumber:\n");
		for (int i = 0; i != halframe; ++i)
		{
			fmt::print("{}  ", count[i]);
			if (count[i] != 0)
			{
				msd[i] /= count[i];
			}
		}
		return msd;
	}
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	/* add some user-defined pargs. */
	real dx = 0.01;
	real dy = 0.01;
	real dz = 0.01;
	real lowx = 0.0;
	real lowy = 0.0;
	real lowz = 0.0;
	real upx = 15.0;
	real upy = 15.0;
	real upz = 20.0;
	real cylCenter1 = 0.0;
	real cylCenter2 = 0.0;
	real dr = 0.01;
	real CR = 0.5;
	real Cr = 0.0;
	int NCM = 0;
	int DIMN = 2;

	hd.pa = {
	{ "-NCM",FALSE, etINT, {&NCM},
	  "center of molecule\n0: mass; 1: geometry; 2: charge"
	},
	{ "-dim",FALSE, etINT, {&DIMN},
	  "dim of bin sets\ndim>2 means getting D in 3-dimension; 0 means yz;1 means xz; 2 means xy;"
	},
	{ "-dx",FALSE, etREAL, {&dx},
	  "the distance of dx(nm)"
	},
	{ "-dy",FALSE, etREAL, {&dy},
	  "the distance of dy(nm)"
	},
	{ "-dz",FALSE, etREAL, {&dz},
	  "the distance of dz(nm)"
	},
	{ "-lowx",FALSE, etREAL, {&lowx},
	  "the region of x(nm)"
	},
	{ "-lowy",FALSE, etREAL, {&lowy},
	  "the region of y(nm)"
	},
	{ "-lowz",FALSE, etREAL, {&lowz},
	  "the region of z(nm)"
	},
	{ "-upx",FALSE, etREAL, {&upx},
	  "the region of x(nm)"
	},
	{ "-upy",FALSE, etREAL, {&upy},
	  "the region of y(nm)"
	},
	{ "-upz",FALSE, etREAL, {&upz},
	  "the region of z(nm)"
	},
	{ "-c1",FALSE, etREAL, {&cylCenter1},
	  "the center of the circle1 (nm)"
	},
	{ "-c2",FALSE, etREAL, {&cylCenter2},
	  "the center of the circle2 (nm)"
	},
	{ "-dr",FALSE, etREAL, {&dr},
	  "the dr bin of R (nm)"
	},
	{ "-CR",FALSE, etREAL, {&CR},
	  "the up region of R (nm)"
	},
	{ "-Cr",FALSE, etREAL, {&Cr},
	  "the low region of R (nm)"
	}
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_msd", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();

	int nbinX, nbinY, nbinZ, nbin1;
	double height, lowPos, upPos;
	if (DIMN == 0)
	{
		height = upx - lowx;
		lowPos = lowx;
		upPos = upx;
	}
	if (DIMN == 1)
	{
		height = upy - lowy;
		lowPos = lowy;
		upPos = upy;
	}
	if (DIMN == 2)
	{
		height = upz - lowz;
		lowPos = lowz;
		upPos = upz;
	}

	nbinX = (upx - lowx) / dx;
	nbinY = (upy - lowy) / dy;
	nbinZ = (upz - lowz) / dz;
	nbin1 = (CR - Cr) / dr;

	itp::matd xcur = hd.initPosc(0);   // the current coordinates 

	std::vector<std::vector<int>> index;
	std::vector<int> tempIndex;
	double tmpPos1, tmpPos2, R;

	std::vector<std::vector<double>> pos;
	std::vector<double> tmpx;

	double tmpTime = 0, dt;
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(xcur, 0);

		tempIndex.clear();
		tmpx.clear();

		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (lowPos <= xcur(i, DIMN) && xcur(i, DIMN) <= upPos)
			{
				if (DIMN == 0)
				{
					tmpPos1 = xcur(i, 1) - cylCenter1;
					tmpPos2 = xcur(i, 2) - cylCenter2;
				}
				if (DIMN == 1)
				{
					tmpPos1 = xcur(i, 0) - cylCenter1;
					tmpPos2 = xcur(i, 2) - cylCenter2;
				}
				if (DIMN == 2)
				{
					tmpPos1 = xcur(i, 0) - cylCenter1;
					tmpPos2 = xcur(i, 1) - cylCenter2;
				}
				R = std::sqrt(tmpPos1 * tmpPos1 + tmpPos2 * tmpPos2);

				if (Cr <= R && R <= CR)
				{
					tempIndex.push_back(i);
					tmpx.push_back(xcur(i, DIMN));
				}
			}
		}

		index.push_back(tempIndex);
		pos.push_back(tmpx);

		dt = hd.fr->time - tmpTime;
		tmpTime = hd.fr->time;
	} while (hd.readNextFrame());

	auto duration = hd.calcDuration(index);

	auto msd = hd.calcMSD(pos, index, duration);

	auto file = hd.openWrite(hd.get_ftp2fn(efXVG));

	for (int i = 0; i != msd.size(); ++i)
	{
		fmt::print(file, "{:8.4f} {:8.4f} \n", i * dt, msd[i]);
	}

	return 0;
}
