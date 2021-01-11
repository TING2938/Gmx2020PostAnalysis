#include <itp/gmx>

void calc_energy(int nmi, int nmj,
	double& vdw, double& cou, const itp::vecd& charge1, const itp::vecd& charge2,
	const itp::boxd& pos1, const itp::boxd& pos2, double Lbox[3], double D, const itp::matd& c6, const itp::matd& c12);

int my_main(int argc, char* argv[])
{
	itp::GmxHandle hd(argc, argv);
	int   DIMN = 2;
	int   nbin = 300;
	float D = 1.0;
	float Cr = 0.0;
	float CR = 0.8;
	float lowPos = 0;
	float upPos = 12;

	hd.pa = {
	  { "-dim",FALSE, etINT, {&DIMN},
		"dim of bin sets\ndim>2 means getting D in 3-dimension,0 X,1 Y,2 Z"
	  },
	  { "-D",FALSE, etREAL, {&D},
		"dielectric constants"
	  },
	  { "-Cr",FALSE, etREAL, {&Cr},
		"Minimum radius of spacing(nm)"
	  },
	  { "-CR",FALSE, etREAL, {&CR},
		"Maximam radius of spacing(nm)"
	  },
	  { "-nbin",FALSE, etINT, {&nbin},
		"number of bin sets for number"
	  },
	  { "-low",FALSE, etREAL, {&lowPos},
		"low position of region of molecule/ion (nm)"
	  },
	  { "-up",FALSE, etREAL, {&upPos},
		"up position of region of molecule/ion (nm)"
	  }
	};

	hd.fnm = {
		{ efXVG, "-o", "", ffOPTWR }
	};

	hd.ngrps = 2;
	hd.init();
	itp::matd c6, c12;
	hd.loadLJParameter(0, 1, c6, c12);
	hd.readFirstFrame();

	/* =================== Main body of code ================== */

	auto pos1 = hd.initPos(0);
	auto pos2 = hd.initPos(1);  // pos of atoms
	auto posc1 = hd.initPosc(0);   // center of mols
	auto posc2 = hd.initPosc(1);

	itp::vecd V(nbin), C(nbin), N(nbin);
	V.fill(0);
	C.fill(0);
	N.fill(0);

	int index = 0;
	int meZ = 0;
	int i, j, k;
	double R, dR;
	double dbin = (upPos - lowPos) / nbin;
	double vdw0, cou0, vdw1, cou1;
	int Number;

	itp::veci nframe(nbin);
	itp::veci numbin(nbin);
	itp::vecd vdw2(nbin);
	itp::vecd cou2(nbin);
	itp::vecd nn2(nbin);
	nframe.fill(0);

	std::string name0 = hd.get_ftp2fn(efXVG);
	if (name0.size() == 4)
	{
		name0 = "Energy_Time_{}_{}.xvg"_format(hd.grpname[0], hd.grpname[1]);
	}
	auto file = hd.openWrite(name0);

	do
	{
		hd.loadPosition(pos1, 0);
		hd.loadPositionCenter(posc1, 0);
		hd.loadPosition(pos2, 1);
		hd.loadPositionCenter(posc2, 1);

		numbin.fill(0);
		vdw2.fill(0);
		cou2.fill(0);
		nn2.fill(0);

		for (i = 0; i != posc1.rows(); ++i)
		{
			if (lowPos <= posc1(i, DIMN) && posc1(i, DIMN) < upPos &&
				1.734 <= posc1(i, 0) && posc1(i, 0) < 3.924 &&
				1.106 <= posc1(i, 1) && posc1(i, 1) < 3.006)
			{
				meZ = int((posc1(i, DIMN) - lowPos) / dbin);
				numbin[meZ]++;

				Number = 0;
				vdw1 = 0;
				cou1 = 0;

				for (j = 0; j != posc2.rows(); ++j)
				{
					R = 0;
					for (k = 0; k != 3; ++k)
					{
						dR = hd.periodicity(posc1(i, k) - posc2(j, k), hd.Lbox[k]);
						R += dR * dR;
					}
					R = sqrt(R);
					if (Cr < R && R < CR)
					{
						calc_energy(i, j, vdw0, cou0, hd.charge[0], hd.charge[1], pos1, pos2, hd.Lbox, D, c6, c12);
						vdw1 += vdw0;
						cou1 += cou0;
						Number++;
					}
				}

				if (Number != 0)
				{
					vdw1 /= Number;
					cou1 /= Number;

					vdw2[meZ] += vdw1;
					cou2[meZ] += cou1;
					nn2[meZ] += Number;
				}
			}
		}

		for (i = 0; i != nbin; ++i)
		{
			if (numbin[i] != 0)
			{
				V[i] += vdw2[i] / numbin[i];
				C[i] += cou2[i] / numbin[i];
				N[i] += nn2[i] / numbin[i];
				nframe[i]++;
			}
		}

		if (numbin[0] != 0)
		{
			fmt::print(file, "{:8.4e} {:8.4e} {:8.4e} {:8.4e} {:8.4e}\n", hd.fr->time, vdw2[0] / numbin[0], cou2[0] / numbin[0], vdw2[0] / numbin[0] + cou2[0] / numbin[0], nn2[0] / numbin[0]);
		}

	} while (hd.readNextFrame());

	for (i = 0; i != nbin; ++i)
	{
		V[i] /= nframe[i];
		C[i] /= nframe[i];
		N[i] /= nframe[i];
	}


	return 0;
}

void calc_energy(int nmi, int nmj,
	double& vdw, double& cou, const itp::vecd& charge1, const itp::vecd& charge2,
	const itp::boxd& pos1, const itp::boxd& pos2, double Lbox[3], double D, const itp::matd& c6, const itp::matd& c12)
{
	double r, dr;
	vdw = 0;
	cou = 0;
	for (int m = 0; m != charge1.size(); ++m)
	{
		for (int n = 0; n != charge2.size(); ++n)
		{
			r = 0;
			for (int k = 0; k != 3; ++k)
			{
				dr = itp::GmxHandle::periodicity(pos1(nmi, m)[k] - pos2(nmj, n)[k], Lbox[k]);
				r += dr * dr;
			}
			r = sqrt(r);
			vdw += (c12(m, n) / pow(r, 12) - c6(m, n) / pow(r, 6));
			cou += 138.9355 / D * charge1[m] * charge2[n] / r;;
		}
	}

}

int main(int argc, char* argv[])
{
	return gmx_run_cmain(argc, argv, &my_main);
}
