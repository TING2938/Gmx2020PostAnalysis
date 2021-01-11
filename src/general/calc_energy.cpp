#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	/**
	 * @brief 计算两个分子之间的相互作用
	 * @param mol1 第一个分子的序号
	 * @param mol2 第二个分子的序号
	 * @param vdw 输出的范德华能
	 * @param cou 输出的静电能
	 * @param D 介电常数
	*/
	void calcEnergyBetweenTwoAtoms(int mol1, int mol2, double& vdw, double& cou, double D);

public:
	itp::matd c6, c12;
	itp::boxd pos1, pos2;
	itp::matd posc1, posc2;
};

gmx_main(calc_energy)
{
	Handle hd(argc, argv);
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

	hd.loadLJParameter(0, 1, hd.c6, hd.c12);
	hd.readFirstFrame();

	/* =================== Main body of code ================== */
	hd.pos1.resize(hd.nmol[0], hd.napm[0]);
	hd.pos2.resize(hd.nmol[1], hd.napm[1]);
	hd.posc1.resize(hd.nmol[0], 3);
	hd.posc2.resize(hd.nmol[1], 3);

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


	do
	{
		hd.loadPosition(hd.pos1, 0);
		hd.loadPositionCenter(hd.posc1, 0);
		hd.loadPosition(hd.pos2, 1);
		hd.loadPositionCenter(hd.posc2, 1);

		numbin.fill(0);
		vdw2.fill(0);
		cou2.fill(0);
		nn2.fill(0);

		for (i = 0; i != hd.nmol[0]; ++i)
		{
			if (lowPos < hd.posc1(i, DIMN) && hd.posc1(i, DIMN) < upPos)
			{
				meZ = int((hd.posc1(i, DIMN) - lowPos) / dbin);
				numbin[meZ]++;

				Number = 0;
				vdw1 = 0;
				cou1 = 0;

				for (j = 0; j != hd.nmol[1]; ++j)
				{
					R = 0;
					for (k = 0; k != 3; ++k)
					{
						dR = hd.periodicity(hd.posc1(i, k) - hd.posc2(j, k), hd.Lbox[k]);
						R += dR * dR;
					}
					R = sqrt(R);
					if (Cr < R && R < CR)
					{
						hd.calcEnergyBetweenTwoAtoms(i, j, vdw0, cou0, D);
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
	} while (hd.readNextFrame());

	for (i = 0; i != nbin; ++i)
	{
		V[i] /= nframe[i];
		C[i] /= nframe[i];
		N[i] /= nframe[i];
	}

	std::string name0 = hd.get_ftp2fn(efXVG);
	if (name0.size() == 4)
	{
		name0 = "Energy-{}-{}-{:.2f}-{:.2f}-{:.2f}-{:.2f}-{}.xvg"_format(hd.grpname[0], hd.grpname[1], Cr, CR, lowPos, upPos, nbin);
	}

	auto file = hd.openWrite(name0);

	for (i = 0; i < nbin; i++)
	{
		fmt::print(file, "{:8.4e} {:8.4e} {:8.4e} {:8.4e} {:8.4e}\n",
			i * dbin, V[i], C[i], V[i] + C[i], N[i]);
	}
	return 0;
}

void Handle::calcEnergyBetweenTwoAtoms(int mol1, int mol2, double& vdw, double& cou, double D)
{
	double r, dr;
	vdw = 0;
	cou = 0;

	for (int m = 0; m != napm[0]; ++m)
	{
		for (int n = 0; n != napm[1]; ++n)
		{
			r = 0;
			for (int k = 0; k != 3; ++k)
			{
				dr = periodicity(pos1(mol1, m)[k] - pos2(mol2, n)[k], Lbox[k]);
				r += dr * dr;
			}
			r = sqrt(r);
			vdw += (c12(m, n) / pow(r, 12) - c6(m, n) / pow(r, 6));
			cou += 138.9355 / D * charge[0][m] * charge[1][n] / r;
		}
	}
}
