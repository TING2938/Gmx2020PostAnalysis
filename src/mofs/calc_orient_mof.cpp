/*  Author     : TING
 *  Date       : 2019/05/05
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate paired and unpair number along z axis in channel simulation.
 */

#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	void getPoint(int type1, int type2, int type3)
	{
		type1--;
		type2--;
		type3--;

		for (int i = 0; i != nmol[0]; ++i)
		{
			for (int m = 0; m != DIM; ++m)
			{
				point1(i, m) = pos(i, type1)[m];
				point2(i, m) = pos(i, type2)[m];
				point3(i, m) = pos(i, type3)[m];
			}
		}
	}

	void calcVector(int nm, double v1[3], double v2[3])
	{
		for (int i = 0; i != 3; ++i)
		{
			v1[i] = itp::periodicity(point1(nm, i) - point2(nm, i), Lbox[i]);
			v2[i] = itp::periodicity(point1(nm, i) - point3(nm, i), Lbox[i]);
		}
	}

	// degree
	double calcAngle(double v1[3], double v2[3])
	{
		double r1 = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
		double r2 = std::sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
		double cosAngle = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (r1 * r2);
		return std::acos(cosAngle) / itp::pi * 180;
	}

	// get the vector by plus production
	void plusProduct(double a[3], double b[3], double c[3])
	{
		c[0] = a[0] + b[0];
		c[1] = a[1] + b[1];
		c[2] = a[2] + b[2];
	}

	// get the vector by cross production
	void crossProduct(double a[3], double b[3], double c[3])
	{
		c[0] = a[1] * b[2] - a[2] * b[1];
		c[1] = a[2] * b[0] - a[0] * b[2];
		c[2] = a[0] * b[1] - a[1] * b[0];
	}

public:
	itp::boxd pos;
	itp::matd posc;
	itp::matd point1, point2, point3, orient;
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
	int nAngle = 90;
	int NCM = 0;
	int DIMN = 2;
	int tp1 = 1;
	int tp2 = 2;
	int tp3 = 3;
	int DIMNOrient = 2;
	int method = 1;

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
	},
	{ "-nA",FALSE, etINT, {&nAngle},
	  "number of bin sets of angle"
	},
		{ "-tp1", FALSE, etINT, {&tp1}, "atom type 1"},
		{ "-tp2", FALSE, etINT, {&tp2}, "atom type 2"},
		{ "-tp3", FALSE, etINT, {&tp3}, "atom type 3"},
		{ "-DIMNOrient",FALSE, etINT, {&DIMNOrient}, "number of bin sets along one dimension" },
		{ "-method",FALSE, etINT, {&method}, "method" }
	};

	hd.fnm = {
		{ efXVG, "-o", "orient", ffWRITE }
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
	Eigen::ArrayXi NCMDens(nbin1);
	NCMDens.fill(0);

	hd.pos.resize(hd.nmol[0], hd.napm[0]);
	hd.posc.resize(hd.nmol[0], 3);
	hd.point1.resize(hd.nmol[0], 3);
	hd.point2.resize(hd.nmol[0], 3);
	hd.point3.resize(hd.nmol[0], 3);

	double tmpPos1, tmpPos2, R;
	int meZ, meA;
	double vec1[3], vec2[3], vecOut[3];
	double wallVec[3] = { 0, 0, 0 };
	double angle, dAngle;
	dAngle = 180 / nAngle;
	wallVec[DIMNOrient] = 1.0;

	hd.orient.resize(nbin1, nAngle);
	hd.orient.fill(0);

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPosition(hd.pos, 0);
		hd.loadPositionCenter(hd.posc, 0, NCM);
		hd.getPoint(tp1, tp2, tp3);

		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (lowPos <= hd.posc(i, DIMN) && hd.posc(i, DIMN) <= upPos)
			{
				if (DIMN == 0)
				{
					tmpPos1 = hd.posc(i, 1) - cylCenter1;
					tmpPos2 = hd.posc(i, 2) - cylCenter2;
				}
				if (DIMN == 1)
				{
					tmpPos1 = hd.posc(i, 0) - cylCenter1;
					tmpPos2 = hd.posc(i, 2) - cylCenter2;
				}
				if (DIMN == 2)
				{
					tmpPos1 = hd.posc(i, 0) - cylCenter1;
					tmpPos2 = hd.posc(i, 1) - cylCenter2;
				}
				R = std::sqrt(tmpPos1 * tmpPos1 + tmpPos2 * tmpPos2);

				if (Cr <= R && R <= CR)
				{
					meZ = (R - Cr) / dr;
					NCMDens[meZ]++;

					hd.calcVector(i, vec1, vec2);

					if (method == 1)
					{
						hd.plusProduct(vec1, vec2, vecOut);
					}
					if (method == 2)
					{
						hd.crossProduct(vec1, vec2, vecOut);
					}

					angle = hd.calcAngle(vecOut, wallVec);
					meA = angle / dAngle;

					hd.orient(meZ, meA)++;

				}
			}
		}
	} while (hd.readNextFrame());

	auto file = hd.openWrite(hd.get_opt2fn("-o"));
	for (int i = 0; i != nAngle; ++i)
	{
		fmt::print(file, "{:8.4e}  ", (i + 0.5) * dAngle);
		for (int j = 0; j != nbin1; ++j)
		{
			if (NCMDens[j] > 0.0)
			{
				hd.orient(j, i) /= NCMDens[j];
			}
			fmt::print(file, "{:8.4e} ", hd.orient(j, i));
		}
		fmt::print(file, "\n");
	}
	return 0;
}
