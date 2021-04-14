/*  Author     : TING
 *  Date       : 2020/01/09
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate msd.
 */

#include <itp/gmx>

using array3d = std::array<real, 3>;

class Handle : public itp::GmxHandleFull
{
public:
	using GmxHandleFull::GmxHandleFull;

	/**
	 * @brief 计算每个离子在区域内能持续过久
	*/
	std::vector<std::vector<int>> calcDuration(const std::vector<std::vector<int>>& allIndex)
	{
		int nframe = allIndex.size();
		int halframe = nframe / 2;
		std::vector<std::vector<int>> dura(halframe);
		std::vector<int> tmpdura;

		for (int i = 0; i != halframe; ++i)
		{
			tmpdura.resize(allIndex[i].size());
			std::fill(tmpdura.begin(), tmpdura.end(), 0);

			for (int j = 0; j != allIndex[i].size(); ++j)
			{
				for (int k = 0; k != halframe; ++k)
				{
					if (itp::contains(allIndex[i + k], allIndex[i][j]))
					{
						tmpdura[j]++;
					}
					else
					{
						break;
					}
				}
			}
			dura[i].swap(tmpdura);
		}
		return dura;
	}

	itp::vecd calcQMSD(const std::vector<std::vector<array3d>>& allPos, const std::vector<int>& calcDim,
		const std::vector<std::vector<int>>& allIndex, const std::vector<std::vector<int>>& duration, int grp)
	{
		int nframe = allPos.size();
		int halframe = nframe / 2;

		itp::matd qmsd(halframe, calcDim.size());
		itp::vecd count(halframe);
		qmsd.fill(0);
		count.fill(0);
		real pos1, pos2;
		int ndx;

		for (int i = 0; i != halframe; ++i)
		{
			for (int j = 0; j != allIndex[i].size(); ++j)
			{
				for (int k = 0; k != duration[i][j]; ++k)
				{
					for (auto&& m : calcDim)
					{
						pos1 = allPos[i][j][m];
						ndx = std::find(allIndex[i+k].begin(), allIndex[i+k].end(), allIndex[i][j]) - allIndex[i+k].begin();
						pos2 = allPos[i + k][ndx][m];
						qmsd(k, m) += totCharge[grp][allIndex[i][j]] * (pos2 - pos1);
					}
					count[k]++;
				}
			}
		}

		itp::vecd outputQMSD(halframe);
		outputQMSD.fill(0);

		for (int i = 0; i != halframe; ++i)
		{
			for (int m = 0; m < 3; m++)
			{
				outputQMSD[i] += std::pow(qmsd(i, m), 2);
			}
			if (count[i] != 0)
				outputQMSD[i] /= count[i];
		}
		return outputQMSD;
	}

	bool bInRegion(int num, const rvec& low, const rvec& up)
	{
		real tmp;
		for (int m = 0; m < 3; m++)
		{
			tmp = posc(num, m);
			while (tmp > Lbox[m])
				tmp -= Lbox[m];
			while (tmp < 0)
				tmp += Lbox[m];
			if (!(low[m] <= tmp && tmp < up[m]))
			{
				return false;
			}
		}
		return true;
	}

public:
	itp::matd posc;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	/* add some user-defined pargs. */
	rvec lowPos = { 0, 0, 0 };
	rvec upPos = { 30, 30, 30 };
	int com = 0;  // center of molecule, 0:mass, 1:geometry, 2:charge
	int qmsdType = 6; // 0:X, 1:Y, 2:Z, 3:XY, 4:YZ, 5:XZ, 6:XYZ 
	std::vector<std::vector<int>> qmsdCalcDim = {
		{XX}, {YY}, {ZZ},
		{XX, YY}, {YY, ZZ}, {XX, ZZ},
		{XX, YY, ZZ},
	};

	hd.pa = {
		{ "-low", FALSE, etRVEC, {lowPos}, "low position" },
		{ "-up", FALSE, etRVEC, {upPos}, "up position" },
		{ "-com", FALSE, etINT, {&com}, "center of molecule, 0:mass, 1:geometry, 2:charge"},
		{ "-type", FALSE, etINT, {&qmsdType}, "type to calculate. 0:X, 1:Y, 2:Z, 3:XY, 4:YZ, 5:XZ, 6:XYZ"},
	};

	hd.fnm = {
		{ efXVG, "-o", "calc_qmsd_Ft", ffWRITE }
	};

	int grp = 0;
	hd.init();
	hd.readFirstFrame();

	hd.initPosc(hd.posc, grp);

	std::vector<std::vector<int>> allIndex;
	std::vector<int> tempIndex;

	std::vector<std::vector<array3d>> allPos;
	std::vector<array3d> tmpPos;

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc, grp, com);
		tempIndex.clear();
		tmpPos.clear();
		// Assume there's no periodicity problem in the pore
		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (hd.bInRegion(i, lowPos, upPos))
			{
				tempIndex.push_back(i);
				tmpPos.push_back({ hd.posc(i, XX), hd.posc(i, YY), hd.posc(i, ZZ) });
			}
		}
		allIndex.push_back(tempIndex);
		allPos.push_back(tmpPos);
	} while (hd.readNextFrame());

	auto duration = hd.calcDuration(allIndex);

	auto msd = hd.calcQMSD(allPos, qmsdCalcDim[qmsdType], allIndex, duration, grp);

	auto file = hd.openWrite(hd.get_ftp2fn(efXVG));

	for (int i = 0; i != msd.size(); ++i)
	{
		fmt::print(file, "{:8.4f} {:8.4f} \n", i * hd.dt, msd[i]);
	}

	return 0;
}
