/*  Author     : TING
 *  Date       : 2019/11/29
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate lifttime of paired number along z axis in channel simulation.
 */

#include <itp/gmx>

struct SparseMatrix
{
	vecu ndx1; // index of group1;
	vecu ndx2; // index of group2;
	vecd dist; // square distance of molecule center;

	void apend(size_t n1, size_t n2, double distance)
	{
		ndx1.append(n1);
		ndx2.append(n2);
		dist.append(distance);
	}

	void clear()
	{
		ndx1.clear();
		ndx2.clear();
		dist.clear();
	}
};

void getDistance(SparseMatrix& dist, matd& pos1, matd& pos2, double Lbox[3], double rSquare);
void getUnique(SparseMatrix& dist, SparseMatrix& unique);

gmx_main(temp)
{
	itp::GmxHandle hd(argc, argv);
	hd.desc = {
		"calculate lifttime of paired number along z axis in channel simulation.\n"
	};

	hd.ngrps = 2; /* grp1 and grp2, anion and cation. */

	/* user-defined pargs. */
	real r1 = 0.5; // nm;
	real r2 = 0.5; // nm;
	real lowPos = 0; // nm;
	real upPos = 30; // nm;

	hd.pa = {
		{ "-r1",   FALSE, etREAL, {&r1},   "radius of mol1" },
		{ "-r2",   FALSE, etREAL, {&r2},   "radius of mol2" },
		{ "-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_pair_lifttime", ffWRITE }
	};

	hd.init();

	real rSquare = (r1 + r2) * (r1 + r2);

	int nm1 = hd.natoms[0].size(); // nr. of molecule
	int nm2 = hd.natoms[1].size();

	hd.readFirstFrame();

	itp::Vector<SparseMatrix> dipole;

	matd pos1(nm1, 3);
	matd pos2(nm2, 3);
	itp::Vector<matd> allPos1, allPos2;

	SparseMatrix distance, unique;

	/* ---------------------------------------------------------------------------- */
	do {
		hd.loadPositionCenter(pos1, 0);
		hd.loadPositionCenter(pos2, 1);
		allPos1.append(pos1);
		allPos2.append(pos2);

		distance.clear();
		unique.clear();
		getDistance(distance, pos1, pos2, hd.Lbox, rSquare);
		getUnique(distance, unique);
		
		dipole.append(unique);

	} while (hd.readNextFrame());

	/* ---------------------------------------------------------------------- */

	int nframe = dipole.size();
	int frame1 = nframe * 0.5;
	int frame2 = nframe * 0.5;
	veci count(frame2, 0);

	for (int i = 0; i != frame1; ++i) // i, frame1
	{
		auto& d1 = dipole[i]; // dipole data in frame i;
		for (int m = 0; m != d1.ndx1.size(); ++m) // m, for each dipole
		{
			double centerPos = (allPos1[i](d1.ndx1[m], ZZ) + allPos2[i](d1.ndx2[m], ZZ)) / 2;

			if (lowPos <= centerPos && centerPos <= upPos)
			{
				for (int j = 0; j != frame2; ++j) // j, frame2
				{

					auto& d2 = dipole[i + j]; // dipole data in frame i+j;

					if (d2.ndx1.contains(d1.ndx1[m]))
					{
						auto n = itp::find(d2.ndx1 == d1.ndx1[m])[0];
						if (d1.ndx2[m] == d2.ndx2[n])
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
	}
	
	auto file = hd.openWrite(hd.get_opt2fn("-o"));
	for (int i = 0; i != count.size(); ++i)
	{
		fmt::print(file, "{:8d} {:8d}\n", i, count[i]);
	}
	
	/*********************************************************************************/
	return 0;
}

void getDistance(SparseMatrix& dist, matd& pos1, matd& pos2, double Lbox[3], double rSquare)
{
	double dR, dRc;
	for (size_t i = 0; i != pos1.nrow(); ++i) {
		for (size_t j = 0; j != pos2.nrow(); ++j) {
			dRc = 0;
			for (size_t k = 0; k != 2; ++k) {
				dR = itp::periodicity(pos1(i, k) - pos2(j, k), Lbox[k]);
				dRc += dR * dR;
			}
			dRc += (pos1(i, ZZ) - pos2(j, ZZ)) * (pos1(i, ZZ) - pos2(j, ZZ));
			if (dRc <= rSquare)
				dist.apend(i, j, dRc);
		}
	}
}

void getUnique(SparseMatrix& dist, SparseMatrix& unique)
{
	constexpr double REMOVE = 10000;
	size_t n, n1, n2;
	while (true) {
		n = dist.dist.argmin();
		n1 = dist.ndx1[n];
		n2 = dist.ndx2[n];
		if (dist.dist[n] > REMOVE-1)
			return;

		unique.apend(n1, n2, dist.dist[n]);

		for (auto&& i : itp::find(dist.ndx1 == n1)) {
			dist.dist[i] = REMOVE;
		}

		for (auto&& i : itp::find(dist.ndx2 == n2)) {
			dist.dist[i] = REMOVE;
		}
	}
}

