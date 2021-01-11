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


	void prepData(itp::matd& xcur, itp::matd& xprev, int grp)
	{
		double hbox[3];
		for (int m = 0; m != DIM; ++m)
		{
			hbox[m] = 0.5 * Lbox[m];
		}

		for (int i = 0; i != nmol[grp]; ++i)
		{
			for (int m = 0; m != DIM; ++m)
			{
				while (xcur(i, m) - xprev(i, m) <= -hbox[m])
				{
					xcur(i, m) += Lbox[m];
				}
				while (xcur(i, m) - xprev(i, m) > hbox[m])
				{
					xcur(i, m) -= Lbox[m];
				}
			}
		}
	}

	itp::vecd calcMSD(const itp::matd& pos)
	{

		int nframe = pos.rows();
		int halframe = nframe / 2;
		int nmol = pos.cols();

		itp::vecd msd = itp::vecd::Zero(halframe);
		itp::veci count = itp::veci::Zero(halframe);
		double dx, dy, dz;

		for (int i = 0; i != halframe; ++i)
		{
			for (int j = 0; j != halframe; ++j)
			{
				for (int k = 0; k != nmol; ++k)
				{

					dx = pos(i + j, k) - pos(j, k);
					// dy = pos[i + j](n[k], YY) - pos[j](n[k], YY);
					// dz = pos[i + j](n[k], ZZ) - pos[j](n[k], ZZ);
					msd[i] += dx * dx;
					count[i]++;
				}
			}
		}

		for (int i = 0; i != halframe; ++i)
		{
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
	itp::matd xprev = hd.initPosc(0);  // the previous coordinates
	hd.loadPositionCenter(xprev, 0);

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

	std::vector<int> indexAll, indexRes;
	indexAll = index[0];
	for (int i = 1; i != index.size(); ++i)
	{
		indexRes.clear();
		std::set_intersection(indexAll.begin(), indexAll.end(), index[i].begin(), index[i].end(),
			std::inserter(indexRes, indexRes.begin()));
		indexAll.swap(indexRes);
	}

	fmt::print("\n\nTotal molecules in pore: {}\n\n", indexAll.size());
	if (!indexAll.empty())
	{
		itp::matd pos1(pos.size(), indexAll.size());

		for (int i = 0; i != pos.size(); ++i)
		{
			for (int j = 0; j != indexAll.size(); ++j)
			{
				auto ndx = std::find(index[i].begin(), index[i].end(), indexAll[j]) - index[i].begin();
				pos1(i, j) = pos[i][ndx];
			}
		}

		auto msd = hd.calcMSD(pos1);

		auto file = hd.openWrite(hd.get_ftp2fn(efXVG));

		for (int i = 0; i != msd.size(); ++i)
		{
			fmt::print(file, "{:8.4f} {:8.4f} \n", i * dt, msd[i]);
		}
	}
	return 0;
}
