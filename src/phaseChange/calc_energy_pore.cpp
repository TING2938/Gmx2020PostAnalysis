#include <itp/gmx>

void calc_energy(int nmi, int nmj,
	double& vdw, double& cou, const vecd& charge1, const vecd& charge2,
	const itp::Vector<matd>& pos1, const itp::Vector<matd>& pos2, double Lbox[3], double D, const matd& c6, const matd& c12);

int my_main(int argc, char *argv[])
{
	itp::GmxHandle hd(argc, argv);
    
    real D = 1.0;
	real Cr = 0.0;
    real CR = 0.8;
	int lateral = 2;
	int dim1 = 0;
	real lowPos1 = -2; // nm;
	real upPos1 = 30; // nm;
	int dim2 = 2;
	real lowPos2 = -2;
	real upPos2 = 30;

	hd.pa = {
	  { "-d1",FALSE, etINT, {&dim1},
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
	  { "-low1",FALSE, etREAL, {&lowPos1}, "low position of region of molecule/ion (nm)"
	  },
	  { "-up1",FALSE, etREAL, {&upPos1}, "up position of region of molecule/ion (nm)" },
	  { "-d2", FALSE, etINT, {&dim2}, "direction to calc. 0(x), 1(y), 2(z)"},
	  { "-up2", FALSE, etREAL, {&upPos2}, "up position of region of molecule/ion (nm)" },
      { "-low2", FALSE, etREAL, {&lowPos2}, "low position of region of molecule/ion (nm)" }
	};


    hd.fnm = {
        { efXVG, "-o", "", ffOPTWR }
    };

	hd.ngrps = 2;
	hd.init();
	matd c6, c12;
	hd.loadLJParameter(0, 1, c6, c12);
	hd.readFirstFrame();
	
    /* =================== Main body of code ================== */

	auto pos1 = hd.initPos(0);
	auto pos2 = hd.initPos(1);  // pos of atoms
	auto posc1 = hd.initPosc(0);   // center of mols
	auto posc2 = hd.initPosc(1);

	vecd V(nbin, 0), C(nbin, 0), N(nbin, 0);

    int index = 0; 
    int meZ = 0;
	int i, j, k;
	double R, dR;
	double vdw0, cou0, vdw1, cou1;
	int Number;

	veci nframe(nbin, 0);
	veci numbin(nbin, 0);
	vecd vdw2(nbin, 0);
	vecd cou2(nbin, 0);
	vecd nn2(nbin, 0);

    do { 
		hd.loadPosition(pos1, 0); 
		hd.loadPositionCenter(posc1, 0);
		hd.loadPosition(pos2, 1);
		hd.loadPositionCenter(posc2, 1);

		numbin.fill(0);
		vdw2.fill(0);
		cou2.fill(0);
		nn2.fill(0);

		for (i = 0; i != hd.natoms[0].size(); ++i) {
			if (lowPos < posc1(i, DIMN) && posc1(i, DIMN) < upPos) {
				meZ = int((posc1(i, DIMN) - lowPos) / dbin);

				Number = 0;
				vdw1 = 0;
				cou1 = 0;

				for (j = 0; j != hd.natoms[1].size(); ++j) {
					R = 0;
					for (k = 0; k != 3; ++k) {
						dR = itp::periodicity(posc1(i, k) - posc2(j, k), hd.Lbox[k]);
						R += dR * dR;
					}
					R = sqrt(R); 
					if (Cr < R && R < CR) { 
						calc_energy(i, j, vdw0, cou0, hd.charge[0][i], hd.charge[1][j], pos1, pos2, hd.Lbox, D, c6, c12);
						vdw1 += vdw0;
						cou1 += cou0;
						Number++;
					}
				} 
				
				if (Number != 0) { 
					vdw1 /= Number;
					cou1 /= Number; 

					vdw2[meZ] += vdw1;
					cou2[meZ] += cou1;
					nn2[meZ] += Number;
					numbin[meZ]++;
				}
			} 
		}

		for (i = 0; i != nbin; ++i) {
			if (numbin[i] != 0) {
				V[i] += vdw2[i] / numbin[i];
				C[i] += cou2[i] / numbin[i];
				N[i] += nn2[i] / numbin[i];
				nframe[i]++;
			}
		}

    } while (hd.readNextFrame());

	for (i = 0; i != nbin; ++i) {
		V[i] /= nframe[i];
		C[i] /= nframe[i];
		N[i] /= nframe[i];
	}

	string name0 = hd.get_ftp2fn(efXVG);
	if (name0.size() == 4) {
		name0 = "Energy2_{}_{}.xvg"_format(hd.grpname[0], hd.grpname[1]);
	}
	
    auto file = hd.openWrite(name0);

    for (i = 0; i < nbin; i++) {
        fmt::print(file, "{:8.4e} {:8.4e} {:8.4e} {:8.4e} {:8.4e}\n", 
			i*dbin, V[i], C[i], V[i] + C[i], N[i]);
    } 
    return 0;
}

void calc_energy(int nmi, int nmj,  
	double& vdw, double& cou, const vecd& charge1, const vecd& charge2,
	const itp::Vector<matd>& pos1, const itp::Vector<matd>& pos2, double Lbox[3], double D, const matd& c6, const matd& c12)
{
	double r, dr;
	vdw = 0;
	cou = 0;
	for (int m = 0; m != charge1.size(); ++m) {
		for (int n = 0; n != charge2.size(); ++n) {
			r = 0;
			for (int k = 0; k != 3; ++k) {
				dr = itp::periodicity(pos1[nmi][m][k] - pos2[nmj][n][k], Lbox[k]);
				r += dr * dr;
			}
			r = sqrt(r);
			vdw += (c12(m, n) / pow(r, 12) - c6(m, n) / pow(r, 6));
			cou += 138.9355 / D * charge1[m] * charge2[n] / r;;
		}
	}

}

int main(int argc, char *argv[])
{
    return gmx_run_cmain(argc, argv, &my_main);
}
