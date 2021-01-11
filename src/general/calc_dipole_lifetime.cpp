/*  Author     : TING
 *  Date       : 2019/11/29
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate lifttime of paired number along z axis in channel simulation.
 */

#include <itp/gmx>

using vecu = std::vector<unsigned int>;
using vecd = std::vector<double>;

struct SparseMatrix
{
	vecu ndx1; // index of group1;
	vecu ndx2; // index of group2;
	vecd dist; // square distance of molecule center;

	void apend(size_t n1, size_t n2, double distance)
	{
		ndx1.push_back(n1);
		ndx2.push_back(n2);
		dist.push_back(distance);
	}

	void clear()
	{
		ndx1.clear();
		ndx2.clear();
		dist.clear();
	}
};

void getDistance(SparseMatrix& dist, itp::matd& pos1, itp::matd& pos2, double Lbox[3], double rSquare);
void getUnique(SparseMatrix& dist, SparseMatrix& unique);

gmx_main(temp)
{
	itp::GmxHandle hd(argc, argv);
	hd.desc = {
		"calculate lifttime of paired number along z axis in channel simulation.\n"
	};

	hd.ngrps = 2; /* grp1 and grp2, anion and cation. */

	/* user-defined pargs. */
	real R = 0.5; // nm
	real lowPos = 0; // nm;
	real upPos = 30; // nm;

	hd.pa = {
		{ "-R",   FALSE, etREAL, {&R},   "R between two types (nm)" },
		{ "-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_lifttime", ffWRITE }
	};

	hd.init();
	
	int nm1 = hd.napm[0]; // nr. of molecule
	int nm2 = hd.napm[1];

	hd.readFirstFrame();

	itp::matd pos1(nm1, 3);
	itp::matd pos2(nm2, 3);
	std::vector<itp::matd> allPos1, allPos2;

	/* ---------------------------------------------------------------------------- */
	do {
		hd.loadPositionCenter(pos1, 0);
		hd.loadPositionCenter(pos2, 1);
		allPos1.push_back(pos1);
		allPos2.push_back(pos2);

	} while (hd.readNextFrame());

	/* ---------------------------------------------------------------------- */

	int nframe = hd.nframe;
	int frame1 = nframe * 0.5;
	int frame2 = nframe * 0.5;

	itp::veci count(frame1);
	count.fill(0);

	for (int ii = 0; ii < frame1; ii++)
	{
		for (int jj = 0; jj < frame2; jj++)
		{
			for (int mm = 0; mm < nm1; mm++)
			{
				for (int nn = 0; mm < nm2; nn++)
				{

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

