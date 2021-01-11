/*  Author     : TING
 *  Date       : 2020/12/04
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate free and bond number in bulk simulation.
 */
#include <itp/gmx>
#include <numeric>

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

		if (distance.value.empty())
		{
            fmt::print("empty occur!\n");
			return;
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
	
	auto calcLifeTime(const std::vector<SparseMatrix>& allInfo)
	{
		int nframe = allInfo.size();
		int halframe = nframe / 2;
		std::vector<int>::const_iterator it;
		std::vector<int> count(halframe);
		std::fill(count.begin(), count.end(), 0);

		for (int i = 0; i != halframe; ++i) // i, frame1
		{
			auto& d1 = allInfo[i]; // dipole data in frame i;
			for (int m = 0; m != d1.row.size(); ++m) // m, for each dipole
			{
				for (int j = 0; j != halframe; ++j) // j, frame2
				{
					auto& d2 = allInfo[i + j]; // dipole data in frame i+j;

					it = std::find(d2.row.begin(), d2.row.end(), d1.row[m]);
					if (it != d2.row.end())
					{
						auto n = it - d2.row.begin();
						if (d1.col[m] == d2.col[n])
						{
							count[j]++;
						}
						else
						{
							break;
						}
					}
					else
					{
						break;
					}
				}
			}
		}
		return count;
	}

public:
	itp::matd posc1;
	itp::matd posc2;
};

gmx_main(calc_bond)
{
	Handle hd(argc, argv);

	hd.ngrps = 2;
	
	/* user-defined pargs. */
	real Rmin = 0.5; // nm;
	int nbin = 100;
	real lowPos = 0.0;
	real upPos = 10.0;

	hd.pa = {
		{ "-Rmin", FALSE, etREAL, {&Rmin},   "Rmin from RDF" },
		{ "-nbin", FALSE, etINT, {&nbin}, "nbin"},
		{ "-low", FALSE, etREAL, {&lowPos}, "low pos"},
		{ "-up", FALSE, etREAL, {&upPos}, "up pos"}
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_pair", ffWRITE }, 
		{ efXVG, "-ob", "bond_pair", ffWRITE},
		{ efXVG, "-lt", "lifetime", ffWRITE}
	};

	hd.init();
	hd.readFirstFrame();

	hd.posc1.resize(hd.nmol[0], 3);
	hd.posc2.resize(hd.nmol[1], 3);
	double dbin = (upPos - lowPos) / nbin;
	itp::vecd bond(nbin);
	itp::vecd free1(nbin);
	itp::vecd free2(nbin);
	bond.fill(0);
	free1.fill(0);
	free2.fill(0);
	
	auto diff = [](const std::vector<int>& vec1, const std::vector<int>& vec2) {
		if (vec2.empty())
		{
			return vec1;
		}
		std::vector<int> vec1_ = vec1;
		std::vector<int> vec2_ = vec2;
		std::sort(vec1_.begin(), vec1_.end());
		std::sort(vec2_.begin(), vec2_.end());
		std::vector<int> ret;
		auto it = std::set_difference(vec1_.begin(), vec1_.end(),
			vec2_.begin(), vec2_.end(),
			std::inserter(ret, ret.begin()));
		return ret;
	};

	const char* testChar = "1234\0 567";

	std::vector<int> allNdx1(hd.nmol[0]);
	std::vector<int> allNdx2(hd.nmol[1]);
	std::iota(allNdx1.begin(), allNdx1.end(), 0);
	std::iota(allNdx2.begin(), allNdx2.end(), 0);

	SparseMatrix pairedMolecule;
	std::vector<SparseMatrix> allPairInfo;
	double p1, p2;
	int meZ;
	double bondNumber = 0.0;
	
	class TestClass {};
	sizeof(TestClass);

	/* ---------------------------------------------------------------------------- */
	do {
		hd.loadPositionCenter(hd.posc1, 0);
		hd.loadPositionCenter(hd.posc2, 1);

		hd.calcBond(pairedMolecule, Rmin);
		allPairInfo.push_back(pairedMolecule);
		
		bondNumber += pairedMolecule.value.size();

		for (int i = 0; i < pairedMolecule.value.size(); i++)
		{
			p1 = hd.posc1(pairedMolecule.row[i], ZZ);
			p2 = hd.posc2(pairedMolecule.col[i], ZZ);
			if (lowPos <= (p1+p2)/2 && (p1+p2)/2 < upPos)
			{
				meZ = ((p1 + p2) / 2 - lowPos) / dbin;
				if (0 <= meZ && meZ < nbin)
					bond[meZ]++;
			}
		}

		for (auto&& i : diff(allNdx1, pairedMolecule.row))
		{
			p1 = hd.posc1(i, ZZ);
			if (lowPos <= p1 && p1 < upPos)
			{
				meZ = (p1 - lowPos) / dbin;
				if (0 <= meZ && meZ < nbin)
					free1[meZ]++;
			}
		}
		
		for (auto&& i : diff(allNdx2, pairedMolecule.col))
		{
			p2 = hd.posc2(i, ZZ);
			if (lowPos <= p2 && p2 < upPos)
			{
				meZ = (p2 - lowPos) / dbin;
				if (0 <= meZ && meZ < nbin)
					free2[meZ]++;
			}
		}
	} while (hd.readNextFrame());

	bondNumber /= hd.nframe;
	fmt::print("\n\nbond: {}, free1: {}, free2: {}\n\n", bondNumber, hd.nmol[0] - bondNumber, hd.nmol[1] - bondNumber);

	auto obfile = hd.openWrite(hd.get_opt2fn("-ob"));
	for (int i = 0; i < pairedMolecule.value.size(); i++)
	{
		fmt::print(obfile, "{} {} {}\n", pairedMolecule.row[i], pairedMolecule.col[i], pairedMolecule.value[i]);
	}

	double dVol = hd.Lbox[XX] * hd.Lbox[YY] * dbin;
	bond /= (hd.nframe * dVol); // (#/nm^3);
	free1 /= (hd.nframe * dVol);
	free2 /= (hd.nframe * dVol);

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));

	fmt::print(ofile, "# grp1: {}, grp2: {}\n", hd.grpname[0], hd.grpname[1]);
	fmt::print(ofile, "# nm1: {}, nm2: {}\n", hd.nmol[0], hd.nmol[1]);
	fmt::print(ofile, "# boxSize: {} x {} x {}\n", hd.Lbox[XX], hd.Lbox[YY], hd.Lbox[ZZ]);
	fmt::print(ofile, "# lowPos: {}, upPos: {}\n", lowPos, upPos);
	fmt::print(ofile, "# Rmin: {}\n", Rmin);
	fmt::print(ofile, "#    z(nm)     bond      free1      free2  (#/nm^3)\n");

	for (int i = 0; i != nbin; ++i) {
		fmt::print(ofile, "{:8.4e} {:8.4e} {:8.4e} {:8.4e}\n", dbin / 2 + dbin * i, bond[i], free1[i], free2[i]);
	}

	return 0;
}
