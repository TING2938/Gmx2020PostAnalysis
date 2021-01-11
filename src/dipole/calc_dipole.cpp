/*  Author     : TING
 *  Date       : 2019/05/05
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate paired and unpair number along z axis in channel simulation.
 */

#include <itp/gmx>

using vecu = std::vector<size_t>;
using vecd = std::vector<double>;
using matd = itp::ArrayXXd;

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

void getDistance(SparseMatrix& dist, matd& pos1, matd& pos2, double Lbox[3], double rSquare);
void getUnique(SparseMatrix& dist, SparseMatrix& unique);

int calc_pair(int argc, char* argv[])
{
	itp::GmxHandle hd(argc, argv);
	hd.desc = {
		"calculate paired and unpair number along z axis in channel simulation.\n"
	};

	hd.ngrps = 2; /* grp1 and grp2, anion and cation. */

	/* user-defined pargs. */
	int nbin = 100;
	real r1 = 0.5; // nm;
	real r2 = 0.5; // nm;
	real lowPos = 0; // nm;
	real upPos = 30; // nm;

	hd.pa = {
		{ "-nbin", FALSE, etINT,  {&nbin}, "nbins."},
		{ "-r1",   FALSE, etREAL, {&r1},   "radius of mol1" },
		{ "-r2",   FALSE, etREAL, {&r2},   "radius of mol2" },
		{ "-upPos", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-lowPos", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_pair", ffWRITE }
	};

	hd.init();

	real rSquare = (r1 + r2) * (r1 + r2);
	
	int nm1 = hd.nmol[0]; // nr. of molecule
	int nm2 = hd.nmol[1];

	printf("\n\n# grp1: %s, grp2: %s\n", hd.grpname[0], hd.grpname[1]);
	printf("# nm1 = %d, nm2 = %d\n", nm1, nm2);
	printf("# lowPos = %f, upPos = %f\n", lowPos, upPos);
	printf("# r1 = %f, r2 = %f\n\n", r1, r2);

	hd.readFirstFrame();

	vecd dipole(nbin, 0);
	vecd free1(nbin, 0);
	vecd free2(nbin, 0);

	matd pos1(nm1, 3);
	matd pos2(nm2, 3);
	SparseMatrix distance, unique;

	double dbin = (upPos - lowPos) / nbin;
	size_t meZ = 0;
	double p1, p2;

	vecu allNdx1 = itp::arange<size_t>(nm1);
	
	vecu allNdx2 = itp::arange<size_t>(nm2);

	/* ---------------------------------------------------------------------------- */
	do {
		hd.loadPositionCenter(pos1, 0);
		hd.loadPositionCenter(pos2, 1);

		distance.clear();
		unique.clear();
		getDistance(distance, pos1, pos2, hd.Lbox, rSquare);
		getUnique(distance, unique);

		for (size_t i = 0; i != unique.dist.size(); ++i) {
			p1 = pos1(unique.ndx1[i], ZZ);
			p2 = pos2(unique.ndx2[i], ZZ);
			if (lowPos <= p1 && p1 <= upPos && lowPos <= p2 && p2 <= upPos) {
				meZ = (size_t)(((p1 + p2) / 2 - lowPos) / dbin);
				++dipole[meZ];
			}
		}

		for (auto&& i : (allNdx1.setDifference(unique.ndx1))) {
			p1 = pos1(i, ZZ);
			if (lowPos <= p1 && p1 <= upPos) {
				meZ = size_t((p1 - lowPos) / dbin);
				++free1[meZ];
			}
		}

		for (auto&& i : (allNdx2.setDifference(unique.ndx2))) {
			p2 = pos2(i, ZZ);
			if (lowPos <= p2 && p2 <= upPos) {
				meZ = (size_t)((p2 - lowPos) / dbin);
				++free2[meZ];
			}
		}

	} while (hd.readNextFrame());

	/* ---------------------------------------------------------------------- */
	double dVol = hd.Lbox[XX] * hd.Lbox[YY] * dbin;
	dipole /= (hd.nframe * dVol); // (#/nm^3);
	free1 /= (hd.nframe * dVol);
	free2 /= (hd.nframe * dVol);

	auto ofile = hd.openWrite(hd.get_opt2fn("-o")); 

	fprintf(ofile, "# grp1: %s, grp2: %s\n", hd.grpname[0], hd.grpname[1]);
	fprintf(ofile, "# nm1: %d, nm2: %d\n", nm1, nm2);
	fprintf(ofile, "# boxSize: %f x %f x %f\n", hd.Lbox[XX], hd.Lbox[YY], hd.Lbox[ZZ]);
	fprintf(ofile, "# lowPos: %f, upPos: %f\n", lowPos, upPos);
	fprintf(ofile, "# r1: %f, r2: %f\n", r1, r2);
	fprintf(ofile, "#    z(nm)     dipole      free1      free2  (#/nm^3)\n"); 

	for (size_t i = 0; i != nbin; ++i) {
		fprintf(ofile, "%8.4e %8.4e %8.4e %8.4e\n", dbin / 2 + dbin * i, dipole[i], free1[i], free2[i]);
	}
	/*********************************************************************************/
	return 0;
}

void getDistance(SparseMatrix& dist, matd& pos1, matd& pos2, double Lbox[3], double rSquare)
{
	double dR, dRc;
	for (size_t i = 0; i != pos1.rows(); ++i) {
		for (size_t j = 0; j != pos2.rows(); ++j) {
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

int main(int argc, char* argv[])
{
	return gmx_run_cmain(argc, argv, &calc_pair);
}
