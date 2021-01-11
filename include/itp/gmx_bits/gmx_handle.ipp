#include "gmx_handle.hpp"


namespace itp
{
	inline GmxHandle::GmxHandle(int argc, char** argv) : argc(argc), argv(argv),
		flags(TRX_READ_X), ngrps(1), nframe(0)
	{
	}

	inline void GmxHandle::init(bool fullMolecule)
	{
		fnm.push_back({ efTRX, "-f", nullptr, ffREAD });
		fnm.push_back({ efTPR, "-s", nullptr, ffREAD });
		fnm.push_back({ efNDX, "-n", nullptr, ffOPTRD });

		if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, (int)fnm.size(),
			fnm.data(), (int)pa.size(), pa.data(), (int)desc.size(), desc.data(), 0, nullptr, &oenv))
		{
			exit(0);
		}

		top = new t_topology;
		ir = new t_inputrec;
		int natoms;
		ePBC = read_tpx_top(get_ftp2fn(efTPR), ir, nullptr, &natoms, nullptr, nullptr, top);

		snew(grpname, ngrps);
		snew(index, ngrps);
		snew(ngx, ngrps);
		snew(fr, 1);

		printf("\n\033[31mSpecify %d group%s to analysis:\33[0m\n", ngrps, (ngrps > 1) ? "s" : "");
		get_index(&top->atoms, get_ftp2fn_null(efNDX), ngrps, ngx, index, grpname);

		if (fullMolecule)
		{
			napm.resize(ngrps);
			nmol.resize(ngrps);
			charge.resize(ngrps);
			mass.resize(ngrps);

			for (int i = 0; i != ngrps; ++i)
			{
				printf("\n\033[31m### Group %d: '%s': \033[0m\n", i, grpname[i]);
				napm[i] = get_natom_per_mol(i);
				nmol[i] = ngx[i] / napm[i];
				printf("# Nr. of atoms in each molecule : %d\n", napm[i]);
				printf("# Nr. of molecules in group: %d\n", nmol[i]);

				mass[i] = get_mass(i);
				charge[i] = get_charge(i);

				printf("# Print infomation of molecule:\n");
				printf("\033[33m# nr.   name     mass   charge\033[0m\n");
				printf("------------------------------\n");
				for (int j = 0; j != napm[i]; ++j)
				{
					printf("%5d%7s%9.3f%9.3f\n", j, *top->atoms.atomname[index[i][j]], mass[i][j], charge[i][j]);
				}
				printf("------------------------------\n");
				printf("%12s%9.3f%9.3f\n", "Total", mass[i].sum(), charge[i].sum());
			}
			printf("\n");
		}
	}

	inline bool GmxHandle::readFirstFrame()
	{
		bool b = read_first_frame(oenv, &status, get_ftp2fn(efTRX), fr, flags);
		b && (nframe = 1);
		Lbox[XX] = fr->box[XX][XX];
		Lbox[YY] = fr->box[YY][YY];
		Lbox[ZZ] = fr->box[ZZ][ZZ];
		time = fr->time;
		preTime = fr->time;
		dt = 0;
		return b;
	}

	inline bool GmxHandle::readNextFrame()
	{
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
		return b;
	}

	inline boxd GmxHandle::initPos(int grp)
	{
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}
		return boxd(nmol[grp], napm[grp]);
	}

	inline matd GmxHandle::initPosc(int grp)
	{
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}
		return matd(nmol[grp], 3);
	}

	inline void GmxHandle::loadPosition(boxd& pos, int grp)
	{
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		size_t i, j, k;
		int molN;
		int n = 0;
		double tmpPos[3], halfLbox[3];

		for (k = 0; k < 3; k++)
			halfLbox[k] = Lbox[k] / 2;

		for (i = 0; i != nmol[grp]; i++)
		{
			molN = index[grp][n];
			for (j = 0; j != napm[grp]; j++)
			{
				for (k = 0; k != 3; k++)
				{
					tmpPos[k] = fr->x[index[grp][n]][k];
					// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
					while (tmpPos[k] - fr->x[molN][k] > halfLbox[k])
						tmpPos[k] -= Lbox[k];
					while (tmpPos[k] - fr->x[molN][k] < -halfLbox[k])
						tmpPos[k] += Lbox[k];
					pos(i, j)[k] = tmpPos[k];
				}
				n++;
			}
		}
	}

	inline void GmxHandle::loadPositionCenter(matd& posc, int grp, int center)
	{
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		int i, j, k, molN;
		int n = 0;
		double tmpPos, tmpCenter, halfLbox[3];

		vecd atomCenter;
		if (center == 0)
		{
			atomCenter = mass[grp];
		}
		else if (center == 1)
		{
			atomCenter.resize(napm[grp]);
			atomCenter.fill(1);
		}
		else if (center == 2)
		{
			atomCenter = charge[grp];
		}

		for (k = 0; k < 3; k++)
			halfLbox[k] = Lbox[k] / 2;

		for (i = 0; i != nmol[grp]; i++)
		{
			molN = index[grp][n];
			for (k = 0; k != 3; k++)
			{
				tmpCenter = 0.0;  // initiate for each molecule
				for (j = 0; j != napm[grp]; j++)
				{
					tmpPos = fr->x[molN + j][k];
					// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
					while (tmpPos - fr->x[molN][k] > halfLbox[k])
						tmpPos -= Lbox[k];
					while (tmpPos - fr->x[molN][k] < -halfLbox[k])
						tmpPos += Lbox[k];
					tmpCenter += tmpPos * atomCenter[j];
				}
				posc(i, k) = tmpCenter / atomCenter.sum();
			}
			n += napm[grp];
		}
	}

	inline void GmxHandle::loadVelocity(boxd& vel, int grp)
	{
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		size_t i, j, k;
		int molN;
		int n = 0;
		double tmpPos[3];

		for (i = 0; i != nmol[grp]; i++)
		{
			molN = index[grp][n];
			for (j = 0; j != napm[grp]; j++)
			{
				for (k = 0; k != 3; k++)
				{
					tmpPos[k] = fr->v[index[grp][n]][k];
					vel(i, j)[k] = tmpPos[k];
				}
				n++;
			}
		}
	}

	inline void GmxHandle::loadVelocityCenter(matd& velc, int grp, int com)
	{
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		int i, j, k, molN;
		int n = 0;
		double tmpPos, tmpCenter;

		vecd atomCenter;
		if (com == 0)
		{
			atomCenter = mass[grp];
		}
		else if (com == 1)
		{
			atomCenter.resize(napm[grp]);
			atomCenter.fill(1);
		}
		else if (com == 2)
		{
			atomCenter = charge[grp];
		}

		for (i = 0; i != nmol[grp]; i++)
		{
			molN = index[grp][n];
			for (k = 0; k != 3; k++)
			{
				tmpCenter = 0.0;  // initiate for each molecule
				for (j = 0; j != napm[grp]; j++)
				{
					tmpPos = fr->v[molN + j][k];
					tmpCenter += tmpPos * atomCenter[j];
				}
				velc(i, k) = tmpCenter / atomCenter.sum();
			}
			n += napm[grp];
		}
	}
	
	inline void GmxHandle::loadForce(boxd& force, int grp)
	{
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		size_t i, j, k;
		int molN;
		int n = 0;
		double tmpPos[3];

		for (i = 0; i != nmol[grp]; i++)
		{
			molN = index[grp][n];
			for (j = 0; j != napm[grp]; j++)
			{
				for (k = 0; k != 3; k++)
				{
					tmpPos[k] = fr->f[index[grp][n]][k];
					force(i, j)[k] = tmpPos[k];
				}
				n++;
			}
		}
	}

	inline void GmxHandle::loadForceCenter(matd& forcec, int grp, int com)
	{
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		int i, j, k, molN;
		int n = 0;
		double tmpPos, tmpCenter;

		vecd atomCenter;
		if (com == 0)
		{
			atomCenter = mass[grp];
		}
		else if (com == 1)
		{
			atomCenter.resize(napm[grp]);
			atomCenter.fill(1);
		}
		else if (com == 2)
		{
			atomCenter = charge[grp];
		}

		for (i = 0; i != nmol[grp]; i++)
		{
			molN = index[grp][n];
			for (k = 0; k != 3; k++)
			{
				tmpCenter = 0.0;  // initiate for each molecule
				for (j = 0; j != napm[grp]; j++)
				{
					tmpPos = fr->f[molN + j][k];
					tmpCenter += tmpPos * atomCenter[j];
				}
				forcec(i, k) = tmpCenter / atomCenter.sum();
			}
			n += napm[grp];
		}
	}

	inline void GmxHandle::loadLJParameter(int group1, int group2, matd& c6, matd& c12)
	{
		if ((group1 >= ngrps) || (group2 >= ngrps))
		{
			gmx_fatal(FARGS, "grp(%d and %d) should less than ngrps(%d)!", group1, group2, ngrps);
		}
		int n1, n2;
		int ntypes = top->atomtypes.nr;
		int type;
		c6.resize(napm[group1], napm[group2]);
		c12.resize(napm[group1], napm[group2]);
		printf("\n\033[31m### LJ Parameters(Group %d and %d):\033[0m\n", group1, group2);
		printf("[  at1  at2]: %12s %12s\n", "c6", "c12");
		printf("---------------------------------------\n");
		int n = 0;
		for (int i = 0; i != napm[group1]; ++i)
		{
			for (int j = 0; j != napm[group2]; ++j)
			{
				n1 = index[group1][i];
				n2 = index[group2][j];
				type = ntypes * (top->atoms.atom[n1].type) + (top->atoms.atom[n2].type);
				c6(i, j) = top->idef.iparams[type].lj.c6;
				c12(i, j) = top->idef.iparams[type].lj.c12;
				if (n++ <= 12)
				{
					printf("[%5s%5s]: %12.4e %12.4e\n", *top->atoms.atomname[n1], *top->atoms.atomname[n2], c6(i, j), c12(i, j));
				}
			}
		}
		if (napm[group1] * napm[group2] > 12)
		{
			printf("(.........)\n");
		}
		printf("\n");
	}

	inline FILE* GmxHandle::openWrite(std::string fnm, bool writeInfo)
	{
		FILE* fp = fopen(fnm.c_str(), "w");
		if (writeInfo)
		{
			char pwd[FILENAME_MAX];
			GetCurrentDir(pwd, sizeof(pwd));
			// get local time
			time_t t;
			char buff[500];
			t = std::time(NULL);
			strftime(buff, 500, "%Y-%m-%d %H:%M:%S %A", localtime(&t));

			fprintf(fp, "# This file was created at %s\n", buff);
			fprintf(fp, "# Working dir : %s\n", pwd);
			fprintf(fp, "# Command line:");
			for (int i = 0; i != argc; ++i)
			{
				fprintf(fp, " %s", argv[i]);
			}
			fprintf(fp, "\n");
		}
		return fp;
	}

	inline const char* GmxHandle::get_ftp2fn(int ftp)
	{
		return ftp2fn(ftp, (int)fnm.size(), fnm.data());
	}

	inline const char* GmxHandle::get_ftp2fn_null(int ftp)
	{
		return ftp2fn_null(ftp, (int)fnm.size(), fnm.data());
	}

	inline const char* GmxHandle::get_opt2fn(const char* opt)
	{
		return opt2fn(opt, (int)fnm.size(), fnm.data());
	}

	inline double GmxHandle::periodicity(double dx, double box)
	{
		while (dx > box / 2.0)
			dx -= box;
		while (dx < -box / 2.0)
			dx += box;
		return dx;
	}

	inline int GmxHandle::get_natom_per_mol(int grp)
	{
		int res, nat, i, j;
		std::vector<int> natom;
		i = 0;
		while (i < ngx[grp])
		{
			res = top->atoms.atom[index[grp][i]].resind;
			nat = 0;
			for (j = top->mols.index[res]; j != top->mols.index[res + 1]; ++j)
			{
				if (i >= ngx[grp] || index[grp][i] != j)
				{
					gmx_fatal(FARGS, "The index group does not consist of whole molecules");
				}
				i++;
				nat++;
			}
			natom.push_back(nat);
		}
		for (auto&& i : natom)
		{
			if (i != natom[0])
			{
				fmt::print(stderr, "\nError: atoms per molecule is not the same in the group!");
				std::exit(-1);
			}
		}
		return natom[0];
	}

	inline vecd GmxHandle::get_mass(int grp)
	{
		vecd mass(napm[grp]);
		for (int i = 0; i != napm[grp]; ++i)
		{
			mass[i] = top->atoms.atom[index[grp][i]].m;
		}
		return mass;
	}

	inline vecd GmxHandle::get_charge(int grp)
	{
		vecd charge(napm[grp]);
		for (int i = 0; i != napm[grp]; ++i)
		{
			charge[i] = top->atoms.atom[index[grp][i]].q;
		}
		return charge;
	}

}
