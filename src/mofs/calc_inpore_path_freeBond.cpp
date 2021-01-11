/*  Author     : TING
 *  Date       : 2019/06/05
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>
#include <set>
#include <vector>
#include <utility>

struct SparseMatrix
{
	std::vector<int> row; // index of row;
	std::vector<int> col; // index of col;
	std::vector<double> value; // value;

	void push_back(int n1, int n2, double val)
	{
		row.push_back(n1);
		col.push_back(n2);
		value.push_back(val);
	}

	void clear()
	{
		row.clear();
		col.clear();
		value.clear();
	}
};

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

	std::vector<int> find(const std::vector<int>& vec, int ndx)
	{
		std::vector<int> ret;
		for (int i = 0; i < vec.size(); i++)
		{
			if (vec[i] == ndx)
				ret.push_back(i);
		}
		return ret;
	}

	void calcBond(SparseMatrix& pair, double Rmin)
	{
		SparseMatrix distance;
		double d = 0.0;
		for (int i = 0; i < nmol[0]; i++)
		{
			for (int j = 0; j < nmol[1]; j++)
			{
				d = calcDistance(i, j);
				if (1e-3 < d && d <= Rmin)
				{
					distance.push_back(i, j, d);
				}
			}
		}
		
		pair.clear();
		constexpr double REMOVE = 100000;
		int n, n1, n2;
		while (true)
		{
			n = std::min_element(distance.value.begin(), distance.value.end()) - distance.value.begin();
			n1 = distance.row[n];
			n2 = distance.col[n];
			if (distance.value[n] > REMOVE - 1)
				return;
			pair.push_back(n1, n2, distance.value[n]);

			for (auto&& i : find(distance.row, n1))
				distance.value[i] = REMOVE;
			for (auto&& i : find(distance.col, n2))
				distance.value[i] = REMOVE;
		}
	}

public:
	itp::matd posc1;
	itp::matd posc2;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	hd.ngrps = 2;

	// add some user-defined pargs.
	real lowPos = 7.761; // nm;
	real upPos = 12.239; // nm;
	real dim = 2;
	real centerX = 2.829;
	real centerY = 2.056;
	real CR = 0.927;
	real Cr = 0.4;
	real Rmin = 0.8; // nm
	real radius = 0.5;

	hd.pa = {
		{ "-upPos", FALSE, etREAL, {&upPos}, "up position of region of pore (nm)" },
		{ "-lowPos", FALSE, etREAL, {&lowPos}, "low position of region of pore (nm)" },
		{ "-dim", FALSE, etREAL, {&dim}, "dim, 0(x), 1(y), 2(z)" },
		{ "-centerX", FALSE, etREAL, {&centerX}, "center X position of mofs pore (nm)" },
		{ "-centerY", FALSE, etREAL, {&centerY}, "center Y position of mofs pore (nm)" },
		{ "-CR", FALSE, etREAL, {&CR}, "Radius of mofs pore (nm)" },
		{ "-Cr", FALSE, etREAL, {&Cr}, "Radius of pore central (nm)" },
		{ "-Rmin", FALSE, etREAL, {&Rmin}, "Rmin in RDF" },
		{ "-radius", FALSE, etREAL, {&radius}, "sum of vdw radius of two types"}
	};

	hd.fnm = {
		{ efXVG, "-o", "inporePath", ffWRITE },
		{ efXVG, "-op", "outputPath", ffWRITE },
		{ efXVG, "-ob", "bond", ffWRITE}
	};
	hd.init();
	hd.readFirstFrame();

	hd.posc1.resize(hd.nmol[0], 3);
	hd.posc2.resize(hd.nmol[1], 3);

	double Rsquare;
	std::vector<std::pair<std::set<int>, std::set<int>>> allPath; // <center, surface>
	std::set<int> tmpCenter, tmpSurface, allIndex;
	
	
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc1, 0);

		tmpCenter.clear();
		tmpSurface.clear();

		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (lowPos <= hd.posc1(i, dim) && hd.posc1(i, dim) < upPos)
			{
				Rsquare = std::pow(hd.posc1(i, XX) - centerX, 2) + std::pow(hd.posc1(i, YY) - centerY, 2);
				if (Rsquare < Cr * Cr)
				{
					tmpCenter.insert(i);
					allIndex.insert(i);
				}
				else if (Rsquare <= CR * CR)
				{
					tmpSurface.insert(i);
					allIndex.insert(i);
				}
			}
		}
		allPath.emplace_back(tmpCenter, tmpSurface);
	} while (hd.readNextFrame());

	fmt::print("tot frame: {}\n", hd.nframe);

	itp::mati pairNumbr(hd.nframe, allIndex.size());
	pairNumbr.fill(0);

	hd.readFirstFrame();
	int frameIndex = 0;
	int molIndex = 0;
	double distance;
	int numB = 0;

	auto fp_path = hd.openWrite(hd.get_opt2fn("-op"));
	auto fp_bond = hd.openWrite(hd.get_opt2fn("-ob"));

	fmt::print(fp_path, "000  ");
	for (auto&& ndx : allIndex)
	{
		fmt::print(fp_path, "{0} {0} ", ndx);
		fmt::print(fp_bond, "{} ", ndx);
	}
	fmt::print(fp_path, "\n");
	fmt::print(fp_bond, "\n");

	SparseMatrix pairedMolecule;
	std::vector<int>::iterator it;

	do
	{
		hd.loadPositionCenter(hd.posc1, 0);
		hd.loadPositionCenter(hd.posc2, 1);

		hd.calcBond(pairedMolecule, radius);

		fmt::print(fp_path, "{:12.8f} ", hd.fr->time);
		fmt::print(fp_bond, "{:12.8f} ", hd.fr->time);

		molIndex = 0;
		for (auto&& ndx : allIndex)
		{
			// for mol path;
			Rsquare = std::pow(hd.posc1(ndx, XX) - centerX, 2) + std::pow(hd.posc1(ndx, YY) - centerY, 2);
			fmt::print(fp_path, "{:10.4f} {:10.4f} ", hd.posc1(ndx, ZZ), std::sqrt(Rsquare));

			// for free and bond
			it = std::find(pairedMolecule.row.begin(), pairedMolecule.row.end(), ndx);
			if (it != pairedMolecule.row.end())
			{
				fmt::print(fp_bond, "{:12.8f} ", pairedMolecule.value[it - pairedMolecule.row.begin()]);
			}
			else
			{
				fmt::print(fp_bond, "{:12.8f} ", 0.0);
			}

			// for pair number;
			numB = 0;
			for (int j = 0; j < hd.nmol[1]; j++)
			{
				distance = hd.calcDistance(ndx, j);
				if (1e-3 <= distance && distance <= Rmin)
				{
					numB++;
				}
			}
			pairNumbr(frameIndex, molIndex) = numB;
			molIndex++;
		}

		fmt::print(fp_path, "\n");
		fmt::print(fp_bond, "\n");
		frameIndex++;
	} while (hd.readNextFrame());


	auto fp = hd.openWrite(hd.get_opt2fn("-o")); // center: 1, surface: 2, outer: 0

	fmt::print(fp, "{:12s} ", "0.000");
	for (auto ndx : allIndex)
	{
		fmt::print(fp, "{:5d} ", ndx);
	}
	for (auto ndx : allIndex)
	{
		fmt::print(fp, "{:5d} ", ndx);
	}
	fmt::print(fp, "\n");

	auto contain = [](const std::set<int>& container, int value) {
		return container.find(value) != container.end();
	};

	fmt::print("tot frame: {}\n", hd.nframe);

	for (int i = 0; i < hd.nframe; i++)
	{
		fmt::print(fp, "{:12.8f} ", hd.dt * i);
		for (auto&& ndx : allIndex)
		{
			if (contain(allPath[i].first, ndx)) // in center
			{
				fmt::print(fp, "{:5d} ", 1);
			}
			else if (contain(allPath[i].second, ndx)) // in surface
			{
				fmt::print(fp, "{:5d} ", 2);
			}
			else // in outer
			{
				fmt::print(fp, "{:5d} ", 0);
			}
		}
		for (int j = 0; j < allIndex.size(); j++)
		{
			fmt::print(fp, "{:5d} ", pairNumbr(i, j));
		}
		fmt::print(fp, "\n");
	}
	
	return 0;
}
