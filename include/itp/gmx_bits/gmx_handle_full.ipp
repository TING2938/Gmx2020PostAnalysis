#include <filesystem>
#include "gmx_handle_full.hpp"

namespace itp
{
	inline GmxHandleFull::GmxHandleFull(int argc, char** argv) : argc(argc), argv(argv),
		flags(TRX_READ_X), ngrps(1), nframe(0), process(0)
	{
		nthreads = omp_get_num_procs();
	}

	inline void GmxHandleFull::init()
	{
		if (process != 0)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		fnm.push_back({ efTRX, "-f", nullptr, ffREAD });
		fnm.push_back({ efTPR, "-s", nullptr, ffREAD });
		fnm.push_back({ efNDX, "-n", nullptr, ffOPTRD });

		pa.push_back({ "-threads", FALSE, etINT, {&nthreads}, "number of threads" });

		if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, (int)fnm.size(),
			fnm.data(), (int)pa.size(), pa.data(), (int)desc.size(), desc.data(), 0, nullptr, &oenv))
		{
			exit(0);
		}
		process = 1;

		top = new t_topology;
		ir = new t_inputrec;
		int natoms;
		ePBC = read_tpx_top(get_ftp2fn(efTPR), ir, nullptr, &natoms, nullptr, nullptr, top);

		snew(grpname, ngrps);
		snew(index, ngrps);
		snew(ngx, ngrps);
		snew(fr, 1);

		fmt::print("\n\033[31mSpecify {} group{} to analysis:\33[0m\n", ngrps, (ngrps > 1) ? "s" : "");
		get_index(&top->atoms, get_ftp2fn_null(efNDX), ngrps, ngx, index, grpname);

		napm.resize(ngrps);
		nmol.resize(ngrps);
		charge.resize(ngrps);
		mass.resize(ngrps);
		totMass.resize(ngrps);
		totCharge.resize(ngrps);

		for (int i = 0; i != ngrps; ++i)
		{
			get_natom_per_mol(i);
			get_mass_charge(i);

			fmt::print("\n\033[31m### Group {}: '{}': \033[0m\n", i, grpname[i]);
			fmt::print("# Nr. of atoms in the first molecule : {}\n", napm[i][0]);
			fmt::print("# Nr. of molecules in group: {}\n", nmol[i]);
			fmt::print("# Print infomation of the first molecule:\n");
			fmt::print("\033[33m# nr.   name     mass   charge\033[0m\n");
			fmt::print("------------------------------\n");
			for (int j = 0; j < napm[i][0]; j++)
			{
				fmt::print("{:>5d}{:>7s}{:>9.3f}{:>9.3f}\n", j, *top->atoms.atomname[index[i][j]], mass[i][0][j], charge[i][0][j]);
			}
			fmt::print("------------------------------\n");
			fmt::print("{:<12s}{:>9.3f}{:>9.3f}\n", "Total", totMass[i][0], totCharge[i][0]);
		}
	}

	inline bool GmxHandleFull::readFirstFrame()
	{
		if (process != 1)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		bool b = read_first_frame(oenv, &status, get_ftp2fn(efTRX), fr, flags);
		b && (nframe = 1);
		Lbox[XX] = fr->box[XX][XX];
		Lbox[YY] = fr->box[YY][YY];
		Lbox[ZZ] = fr->box[ZZ][ZZ];
		time = fr->time;
		preTime = fr->time;
		dt = 0;
		process = 2;
		return b;
	}

	inline bool GmxHandleFull::readNextFrame()
	{
		if (process != 2)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		bool b = read_next_frame(oenv, status, fr);
		b && (++nframe);
		Lbox[XX] = fr->box[XX][XX];
		Lbox[YY] = fr->box[YY][YY];
		Lbox[ZZ] = fr->box[ZZ][ZZ];
		if (b)
		{
			dt = time - preTime;
			preTime = time;
			time = fr->time;
		}
		process = 3;
		return b;
	}

	inline void GmxHandleFull::initPos(boxd& pos, int grp)
	{
		if (process < 1)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}
		int maxNapm = *std::max_element(napm[grp].begin(), napm[grp].end());
		pos.resize(nmol[grp], maxNapm);
	}

	inline void GmxHandleFull::initPosc(matd& posc, int grp)
	{
		if (process < 1)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}
		posc.resize(nmol[grp], 3);
		posc.fill(0);
	}

	inline void GmxHandleFull::loadPosition(boxd& pos, int grp)
	{
		if (process < 2)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		size_t i, j, k;
		int molN;
		int n = 0;
		double tmpPos, halfLbox[3];

		for (k = 0; k < 3; k++)
			halfLbox[k] = Lbox[k] / 2;

		for (i = 0; i < nmol[grp]; i++)
		{
			molN = index[grp][n];
			for (j = 0; j < napm[grp][i]; j++)
			{
				for (k = 0; k < 3; k++)
				{
					tmpPos = fr->x[index[grp][n]][k];
					// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
					while (tmpPos - fr->x[molN][k] > halfLbox[k])
						tmpPos -= Lbox[k];
					while (tmpPos - fr->x[molN][k] < -halfLbox[k])
						tmpPos += Lbox[k];
					pos(i, j)[k] = tmpPos;
				}
				n++;
			}
		}
	}

	inline void GmxHandleFull::loadPositionCenter(matd& posc, int grp, int center)
	{
		if (process < 2)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		int i, j, k, molN;
		int n = 0;
		double tmpPos, tmpCenter, halfLbox[3];

		std::vector<std::vector<double>> atomCenter;
		std::vector<double> centerSum;

		if (center == 0)
		{
			atomCenter = mass[grp];
			centerSum = totMass[grp];
		}
		else if (center == 1)
		{
			centerSum.resize(nmol[grp]);
			for (int i = 0; i < nmol[grp]; i++)
			{
				atomCenter[i].resize(napm[grp][i], 1);
				centerSum[i] = napm[grp][i];
			}
		}
		else if (center == 2)
		{
			atomCenter = charge[grp];
			centerSum = totCharge[grp];
		}

		for (k = 0; k < 3; k++)
			halfLbox[k] = Lbox[k] / 2;

		for (i = 0; i != nmol[grp]; i++)
		{
			molN = index[grp][n];
			for (k = 0; k != 3; k++)
			{
				tmpCenter = 0.0;  // initiate for each molecule
				for (j = 0; j != napm[grp][i]; j++)
				{
					tmpPos = fr->x[molN + j][k];
					// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
					while (tmpPos - fr->x[molN][k] > halfLbox[k])
						tmpPos -= Lbox[k];
					while (tmpPos - fr->x[molN][k] < -halfLbox[k])
						tmpPos += Lbox[k];
					tmpCenter += tmpPos * atomCenter[i][j];
				}
				posc(i, k) = tmpCenter / centerSum[i];
			}
			n += napm[grp][i];
		}
	}

	inline void GmxHandleFull::loadVelocityCenter(matd& velc, int grp, int com)
	{
		if (process < 2)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		int i, j, k, molN;
		int n = 0;
		double tmpPos, tmpCenter;

		std::vector<std::vector<double>> atomCenter;
		std::vector<double> centerSum;
		if (com == 0)
		{
			atomCenter = mass[grp];
			centerSum = totMass[grp];
		}
		else if (com == 1)
		{
			centerSum.resize(nmol[grp]);
			for (int i = 0; i < nmol[grp]; i++)
			{
				atomCenter[i].resize(napm[grp][i], 1);
				centerSum[i] = napm[grp][i];
			}
		}
		else if (com == 2)
		{
			atomCenter = charge[grp];
			centerSum = totCharge[grp];
		}

		for (i = 0; i != nmol[grp]; i++)
		{
			molN = index[grp][n];
			for (k = 0; k != 3; k++)
			{
				tmpCenter = 0.0;  // initiate for each molecule
				for (j = 0; j != napm[grp][i]; j++)
				{
					tmpPos = fr->v[molN + j][k];
					tmpCenter += tmpPos * atomCenter[i][j];
				}
				velc(i, k) = tmpCenter / centerSum[i];
			}
			n += napm[grp][i];
		}
	}

	inline FILE* GmxHandleFull::openWrite(std::string fnm, bool writeInfo)
	{
		if (process < 1)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		FILE* fp = fopen(fnm.c_str(), "w");
		if (writeInfo)
		{		
			fmt::print(fp, "# This file was created at {}\n", itp::localTime());
			fmt::print(fp, "# Working dir : {}\n", std::filesystem::current_path());
			fmt::print(fp, "# Command line:");
			for (int i = 0; i != argc; ++i)
			{
				fmt::print(fp, " {}", argv[i]);
			}
			fmt::print(fp, "\n");
		}
		return fp;
	}

	inline const char* GmxHandleFull::get_ftp2fn(int ftp)
	{
		if (process < 1)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		return ftp2fn(ftp, (int)fnm.size(), fnm.data());
	}

	inline const char* GmxHandleFull::get_ftp2fn_null(int ftp)
	{
		if (process < 1)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		return ftp2fn_null(ftp, (int)fnm.size(), fnm.data());
	}

	inline const char* GmxHandleFull::get_opt2fn(const char* opt)
	{
		if (process < 1)
		{
			gmx_fatal(FARGS, "The function `%s` was called incorrectly!", __FUNCTION__);
		}
		return opt2fn(opt, (int)fnm.size(), fnm.data());
	}

	inline double GmxHandleFull::periodicity(double dx, double box)
	{
		while (dx > box / 2.0)
			dx -= box;
		while (dx < -box / 2.0)
			dx += box;
		return dx;
	}

	inline void GmxHandleFull::get_natom_per_mol(int grp)
	{
		int firstResi = top->atoms.atom[index[grp][0]].resind;
		int j = 0;
		napm[grp].clear();
		napm[grp].push_back(0);

		for (int i = 0; i < ngx[grp]; i++)
		{
			if (top->atoms.atom[index[grp][i]].resind == firstResi)
			{
				napm[grp][j]++;
			}
			else
			{
				firstResi = top->atoms.atom[index[grp][i]].resind;
				j++;
				napm[grp].push_back(1);
			}
		}
		nmol[grp] = j + 1;
	}

	inline void GmxHandleFull::get_mass_charge(int grp)
	{
		int ndx = 0;
		mass[grp].resize(nmol[grp]);
		charge[grp].resize(nmol[grp]);
		totMass[grp].resize(nmol[grp], 0.0);
		totCharge[grp].resize(nmol[grp], 0.0);

		for (int i = 0; i < nmol[grp]; i++)
		{
			mass[grp][i].resize(napm[grp][i], 0.0);
			charge[grp][i].resize(napm[grp][i], 0.0);
			for (int j = 0; j < napm[grp][i]; j++)
			{
				mass[grp][i][j] = top->atoms.atom[index[grp][ndx]].m;
				charge[grp][i][j] = top->atoms.atom[index[grp][ndx]].q;
				totMass[grp][i] += mass[grp][i][j];
				totCharge[grp][i] += charge[grp][i][j];
				ndx++;
			}
		}
	}

}
