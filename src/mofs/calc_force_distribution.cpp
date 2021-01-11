/*  Author     : TING
 *  Date       : 2019/06/05
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>
#include <set>
#include <vector>
#include <utility>
#include <optional>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	std::optional<int> getPbcMe(int mol, rvec lowPos, rvec upPos, rvec nbin)
	{
		std::vector<int> me;
		int              meIndex[3];
		real             postmp[3];
		for (int i = 0; i < 3; i++)
		{
			postmp[i] = posc(mol, i);
		}
		for (int i = 0; i < 3; i++)
		{
			if (postmp[i] + Lbox[i] < upPos[i])
				postmp[i] += Lbox[i];
			if (postmp[i] - Lbox[i] >= lowPos[i])
				postmp[i] -= Lbox[i];
			if (lowPos[i] <= postmp[i] && postmp[i] < upPos[i])
				meIndex[i] = (postmp[i] - lowPos[i]) / (upPos[i] - lowPos[i]) * nbin[i];
			else
			{
				return std::nullopt;
			}
		}
		return meIndex[XX] + meIndex[YY] * nbin[XX] + meIndex[ZZ] * nbin[XX] * nbin[YY];
	}

	void loadForceCenter(itp::matd& velc, int grp, int com=0)
	{
		if (grp >= ngrps)
		{
			gmx_fatal(FARGS, "grp(%d) should less than ngrps(%d)!", grp, ngrps);
		}

		int i, j, k, molN;
		int n = 0;
		double tmpPos, tmpCenter;

		itp::vecd atomCenter;
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
				velc(i, k) = tmpCenter / atomCenter.sum();
			}
			n += napm[grp];
		}
	}

public:
	itp::matd posc;
	itp::matd velc;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	hd.flags = TRX_NEED_X | TRX_NEED_F;
	hd.ngrps = 1;

	// add some user-defined pargs.
	rvec nbin = { 0, 0, 0 };
	rvec lowPos = { 0, 0, 0 };
	rvec upPos = { 100, 100, 100 };
	int ncm = 0;

	hd.pa = {
		{ "-nbin", FALSE, etRVEC, {&nbin}, "nbin" },
		{ "-low", FALSE, etRVEC, {&lowPos}, "low position of region of pore (nm)" },
		{ "-up", FALSE, etRVEC, {&upPos}, "up position of region of pore (nm)" },
		{ "-ncm", FALSE, etINT, {&ncm}, "which center to use, 0(mass), 1(geo), 2(charge)"}
	};

	hd.fnm = {
		{ efXVG, "-o", "force3d", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();
	
	printf("\n");
	rvec Lbox;
	for (int i = 0; i < 3; i++)
	{
		Lbox[i] = hd.fr->box[i][i];
		if (lowPos[i] < 0)
			lowPos[i] = 0;
		if (upPos[i] > Lbox[i] || upPos[i] <= lowPos[i])
			upPos[i] = Lbox[i];
		printf("lowPos[%d]: %.3f, upPos[%d]: %.3f, Lbox[%d]: %.3f\n", i, lowPos[i], i, upPos[i], i, Lbox[i]);
	}
	
	hd.posc.resize(hd.nmol[0], 3);
	hd.velc.resize(hd.nmol[0], 3);
	
	itp::matd force(int(nbin[0] * nbin[1] * nbin[2]), 3);
	force.fill(0);

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc, 0, ncm);
		hd.loadForceCenter(hd.velc, 0, ncm);

		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (auto&& me = hd.getPbcMe(i, lowPos, upPos, nbin); me.has_value())
			{
				for (int k = 0; k < 3; k++)
				{
					force(me.value(), k) += hd.velc(i, k);
				}
			}
		}
		
	} while (hd.readNextFrame());

	fmt::print("tot frame: {}\n", hd.nframe);

	force /= hd.nframe;
	
	auto ofile = hd.openWrite(hd.get_opt2fn("-o"));
	for (int i = 0; i < force.rows(); i++)
	{
		fmt::print(ofile, "{:15.8e} {:15.8e} {:15.8e}\n", force(i, 0), force(i, 1), force(i, 2));
	}
	return 0;
}
