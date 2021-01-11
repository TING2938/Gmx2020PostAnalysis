/*  Author     : TING
 *  Date       : 2019/05/05
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate paired and unpair number along z axis in channel simulation.
 */

#include <itp/gmx>
#include <omp.h>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	double safe_acos(double value)
	{
		if (value <= -1.0)
		{
			return itp::pi;
		}
		else if (value >= 1.0)
		{
			return 0;
		}
		else
		{
			return std::acos(value);
		}
	}

	double calcVolume(double RMOF, double R, double L);

	double calcDistance(int i, int j)
	{
		double ret = 0.0;
		for (int k = 0; k < 3; k++)
		{
			ret += std::pow(posc1(i, k) - posc2(j, k), 2);
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
	int dim = 2; // 0(x), 1(y), 2(z);
	real lowPos = 7.073; // nm;
	real upPos = 12.927; // nm;
	real centerX = 2.829; // nm
	real centerY = 2.056; // nm
	real RMOF = 0.655; // nm
	real RMOF_max = 0.8; // nm
	int NCM = 0; // kind of center
	real dR = 0.002; // nm

	hd.pa = {
		{ "-dim", FALSE, etINT, {&dim}, "dim, 0(x), 1(y), 2(z)"},
		{ "-up", FALSE, etREAL, {&upPos}, "up position of region of molecule/ion (nm)" },
		{ "-low", FALSE, etREAL, {&lowPos}, "low position of region of molecule/ion (nm)" },
		{ "-centerX", FALSE, etREAL, {&centerX}, "center X (nm)"},
		{ "-centerY", FALSE, etREAL, {&centerY}, "center Y (nm)"},
		{ "-Rmof", FALSE, etREAL, {&RMOF}, "radius of MOF pore (nm)"},
		{ "-RmofMax", FALSE, etREAL, {&RMOF_max}, "max radius of MOF pore (nm)"},
		{ "-ncm", FALSE, etINT, {&NCM}, "kind of center, 0(mass), 1(geometry), 2(charge)"},
		{ "-dR", FALSE, etREAL, {&dR}, "dR (nm)"}
	};

	hd.fnm = {
		{ efXVG, "-o", "rdf_pore", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();
	hd.posc1.resize(hd.nmol[0], 3);
	hd.posc2.resize(hd.nmol[1], 3);

	double LMOF = upPos - lowPos;
	double Rmax = LMOF / 2;
	
	int Rbin = Rmax / dR;
	
	itp::vecd rdf(Rbin);
	itp::veci rdfCount(Rbin);
	rdf.fill(0);
	rdfCount.fill(0);
	double L1, L2, rmax, distance;
	int rbin, meR, meL;
	itp::vecd dens;
	itp::vecd volume;
	volume.resize(Rbin);
	volume.fill(0);

	double dL = 0.001;
	int Lbin = RMOF / dL;
	itp::matd totVolume(Lbin, Rbin);
	totVolume.fill(0);

	fmt::print("\ncalculate volume ...\n");

#pragma omp parallel for 
	for (int i = 0; i < Lbin; i++)
	{
		for (int j = 0; j < Rbin; j++)
		{
			totVolume(i, j) = hd.calcVolume(RMOF, (j+1) * dR, (i+1) * dL);
		}
	}

	for (int i = 0; i < Lbin; i++)
	{
		for (int j = Rbin-1; j >0; j--)
		{
			totVolume(i, j) -= totVolume(i, j - 1);
		}
	}
	fmt::print("calculate volume finished ...\n");

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc1, 0);
		hd.loadPositionCenter(hd.posc2, 1);

		for (int i = 0; i < hd.nmol[0]; i++)
		{
			if (lowPos <= hd.posc1(i, dim) && hd.posc1(i, dim) <= upPos)
			{
				L1 = std::sqrt(std::pow(hd.posc1(i, XX) - centerX, 2) + std::pow(hd.posc1(i, YY) - centerY, 2));
				if (L1 < RMOF_max)
				{
					rmax = std::min(hd.posc1(i, dim) - lowPos, upPos - hd.posc1(i, dim));
					if (rmax > 0)
					{
						rbin = rmax / dR;
						dens.resize(rbin);
						dens.fill(0);

						for (int j = 0; j < hd.nmol[1]; j++)
						{
							if (lowPos <= hd.posc2(j, dim) && hd.posc2(j, dim) <= upPos)
							{
								L2 = std::sqrt(std::pow(hd.posc2(j, XX) - centerX, 2) - std::pow(hd.posc2(j, YY) - centerY, 2));
								if (L2 < RMOF_max)
								{
									distance = hd.calcDistance(i, j);
									if (distance <= rmax)
									{
										meR = distance / rmax * rbin;
										dens[meR]++;
									}
								}
							}
						}
						meL = L1 / RMOF * Lbin;
						if (meL >= Lbin)
							meL = Lbin - 1;
						for (int m = 0; m < rbin; m++)
						{
							dens[m] /= totVolume(meL, m);
							rdf[m] += dens[m];
							rdfCount[m]++;
						}
					}
				}
			}
		}
	} while (hd.readNextFrame());

	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	for (int m = 0; m < Rbin; m++)
	{
		if (rdfCount[m] == 0)
			rdf[m] = 0;
		else
			rdf[m] /= rdfCount[m];

		fmt::print(ofile, "{:12.5f}\t{:12.5f}\n", m * dR, rdf[m]);
	}
	return 0;
}

double Handle::calcVolume(double RMOF, double R, double L)
{
	double ret = 0;
	if (R <= RMOF - L)
	{
		ret = 4 * itp::pi * R * R * R / 3;
	}
	else if (R > RMOF + L)
	{
		double hmax = std::sqrt(R * R - (RMOF - L) * (RMOF - L));
		double hmin = std::sqrt(R * R - (RMOF + L) * (RMOF + L));
		double dh = 1e-4;
		double r, alpha, belta, cosAlpha, cosBelta;
		for (double h = hmin; h <= hmax; h += dh)
		{
			r = std::sqrt(R * R - h * h);
			cosAlpha = (r * r + L * L - RMOF * RMOF) / (2.0 * r * L);
			cosBelta = (L * L + RMOF * RMOF - r * r) / (2.0 * L * RMOF);
			alpha = safe_acos(cosAlpha);
			belta = safe_acos(cosBelta);
			ret += dh * (alpha * r * r + belta * RMOF * RMOF - L * RMOF * std::sin(belta));
		}
		ret += itp::pi / 3 * (3 * R - (R - hmax)) * (R - hmax) * (R - hmax);
		ret += hmin * itp::pi * RMOF * RMOF;
		ret *= 2;
	}
	else
	{
		double hmax = std::sqrt(R * R - (RMOF - L) * (RMOF - L));
		double dh = 1e-4;
		double r, alpha, belta, cosAlpha, cosBelta;
		for (double h = 0; h <= hmax; h += dh)
		{
			r = std::sqrt(R * R - h * h);
			cosAlpha = (r * r + L * L - RMOF * RMOF) / (2.0 * r * L);
			cosBelta = (L * L + RMOF * RMOF - r * r) / (2.0 * L * RMOF);
			alpha = safe_acos(cosAlpha);
			belta = safe_acos(cosBelta);
			ret += dh * (alpha * r * r + belta * RMOF * RMOF - L * RMOF * std::sin(belta));
		}
		ret += itp::pi / 3 * (3 * R - (R - hmax)) * (R - hmax) * (R - hmax);
		ret *= 2;
	}
	return ret;
}

