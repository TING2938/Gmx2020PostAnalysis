/*  Author     : TING
 *  Date       : 2019/05/05
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate paired and unpair number along z axis in channel simulation.
 */

#include <itp/gmx>
#include <map>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	double calcDistance(int i, int j)
	{
		double ret = 0.0;
		for (int k = 0; k < 3; k++)
		{
			ret += std::pow(periodicity(posc1(i, k) - posc2(j, k), Lbox[k]), 2);
		}
		return std::sqrt(ret);
	}

public:
	itp::matd posc1;
	itp::matd posc2;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	hd.ngrps = 2; /* number of group(s), grp1 and grp2, anion and cation. */

	/* add some user-defined pargs. */
	int NCM = 0; // kind of center
	real Rmin = 0.8; // nm
	int dim = 2;
	real lowPos = 0.0;
	real upPos = 12.0;

	hd.pa = {
		{ "-ncm", FALSE, etINT, {&NCM}, "kind of center, 0(mass), 1(geometry), 2(charge)"},
		{ "-Rmin", FALSE, etREAL, {&Rmin}, "the distance for which the RDF shows a first minimum (nm)"},
		{ "-dim", FALSE, etINT, {&dim}, "dim to analysis, 0(x), 1(y), 2(z)"},
		{ "-low", FALSE, etREAL, {&lowPos}, "low position (nm)"},
		{ "-up", FALSE, etREAL, {&upPos}, "up position (nm)"}
	};

	hd.fnm = {
		{ efXVG, "-o", "pair_dist", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();
	hd.posc1.resize(hd.nmol[0], 3);
	hd.posc2.resize(hd.nmol[1], 3);

	double L1, L2, distance;
	int pairNumber, numA = 0;
	std::map<int, double> totPair, tmpPair;
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc1, 0, NCM);
		hd.loadPositionCenter(hd.posc2, 1, NCM);

		numA = 0;
		tmpPair.clear();
		for (int i = 0; i < hd.nmol[0]; i++)
		{
			if (lowPos <= hd.posc1(i, dim) && hd.posc1(i, dim) <= upPos)
			{
				numA++;
				pairNumber = 0;
				for (int j = 0; j < hd.nmol[1]; j++)
				{
					distance = hd.calcDistance(i, j);
					if (distance > 1e-3 && distance <= Rmin)
					{
						pairNumber++;
					}
				}
				tmpPair[pairNumber]++;
			}
		}

		for (auto p : tmpPair)
		{
			totPair[p.first] += p.second / numA;
		}

	} while (hd.readNextFrame());

	double totol = 0;
	for (auto&& p : totPair)
	{
		totol += p.second;
	}

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	for (auto&& p : totPair)
	{
		fmt::print(ofile, "{:8d}\t{:12.5f}\n", p.first, p.second / totol);
	}
	return 0;
}



