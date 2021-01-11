/*  Author     : TING
 *  Date       : 2019/06/05
 *  Email      : yeting2938@hust.edu.cn
 */

#include <itp/gmx>
#include <set>
#include <vector>
#include <utility>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	double calcDistance(int i, int j)
	{
		double ret = 0.0;
		for (int k = 0; k < 3; k++)
		{
			ret += std::pow(periodicity(posc1(i, k) - posc2(j, k), Lbox[k]), 2);
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

	hd.ngrps = 2;

	// add some user-defined pargs.
	real lowPos = 7.761; // nm;
	real upPos = 12.239; // nm;
	real dim = 2;
	real centerX = 2.829;
	real centerY = 2.056;
	real CR = 0.927;
	real Cr = 0.4;
	real Rmin = 0.8; // nm

	hd.pa = {
		{ "-upPos", FALSE, etREAL, {&upPos}, "up position of region of pore (nm)" },
		{ "-lowPos", FALSE, etREAL, {&lowPos}, "low position of region of pore (nm)" },
		{ "-dim", FALSE, etREAL, {&dim}, "dim, 0(x), 1(y), 2(z)" },
		{ "-centerX", FALSE, etREAL, {&centerX}, "center X position of mofs pore (nm)" },
		{ "-centerY", FALSE, etREAL, {&centerY}, "center Y position of mofs pore (nm)" },
		{ "-CR", FALSE, etREAL, {&CR}, "Radius of mofs pore (nm)" },
		{ "-Cr", FALSE, etREAL, {&Cr}, "Radius of pore central (nm)" },
		{ "-Rmin", FALSE, etREAL, {&Rmin}, "Rmin in RDF" }
	};

	hd.fnm = {
		{ efXVG, "-o", "inporePath", ffWRITE },
		{ efXVG, "-op", "outputPath", ffWRITE }
	};
	hd.init();
	hd.readFirstFrame();

	hd.posc1.resize(hd.nmol[0], 3);
	hd.posc2.resize(hd.nmol[1], 3);

	double Rsquare;
	std::vector<std::pair<std::set<int>, std::set<int>>> allPath; // <center, surface>
	std::set<int> tmpCenter, tmpSurface, allIndex;
	
	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadPositionCenter(hd.posc1, 0);

		tmpCenter.clear();
		tmpSurface.clear();

		for (int i = 0; i != hd.nmol[0]; ++i)
		{
			if (lowPos <= hd.posc1(i, dim) && hd.posc1(i, dim) < upPos)
			{
				Rsquare = std::pow(hd.posc1(i, XX) - centerX, 2) + std::pow(hd.posc1(i, YY) - centerY, 2);
				if (Rsquare < Cr * Cr)
				{
					tmpCenter.insert(i);
					allIndex.insert(i);
				}
				else if (Rsquare <= CR * CR)
				{
					tmpSurface.insert(i);
					allIndex.insert(i);
				}
			}
		}
		allPath.emplace_back(tmpCenter, tmpSurface);
	} while (hd.readNextFrame());

	fmt::print("tot frame: {}\n", hd.nframe);

	itp::mati pairNumbr(hd.nframe, allIndex.size());
	pairNumbr.fill(0);

	hd.readFirstFrame();
	int frameIndex = 0;
	int molIndex = 0;
	double distance;
	int numB = 0;

	auto fp_path = hd.openWrite(hd.get_opt2fn("-op"));
	fmt::print(fp_path, "000  ");
	for (auto&& ndx : allIndex)
	{
		fmt::print(fp_path, "{0} {0} ", ndx);
	}
	fmt::print(fp_path, "\n");

	do
	{
		hd.loadPositionCenter(hd.posc1, 0);
		hd.loadPositionCenter(hd.posc2, 1);

		fmt::print(fp_path, "{:12.8f} ", hd.fr->time);

		molIndex = 0;
		for (auto ndx : allIndex)
		{
			Rsquare = std::pow(hd.posc1(ndx, XX) - centerX, 2) + std::pow(hd.posc1(ndx, YY) - centerY, 2);
			fmt::print(fp_path, "{:10.4f} {:10.4f} ", hd.posc1(ndx, ZZ), std::sqrt(Rsquare));

			numB = 0;
			for (int j = 0; j < hd.nmol[1]; j++)
			{
				distance = hd.calcDistance(ndx, j);
				if (1e-3 <= distance && distance <= Rmin)
				{
					numB++;
				}
			}
			pairNumbr(frameIndex, molIndex) = numB;
			molIndex++;
		}

		fmt::print(fp_path, "\n");
		frameIndex++;
	} while (hd.readNextFrame());


	auto fp = hd.openWrite(hd.get_opt2fn("-o")); // center: 1, surface: 2, outer: 0

	fmt::print(fp, "{:12s} ", "0.000");
	for (auto ndx : allIndex)
	{
		fmt::print(fp, "{:5d} ", ndx);
	}
	for (auto ndx : allIndex)
	{
		fmt::print(fp, "{:5d} ", ndx);
	}
	fmt::print(fp, "\n");

	auto contain = [](const std::set<int>& container, int value) {
		return container.find(value) != container.end();
	};

	fmt::print("tot frame: {}\n", hd.nframe);

	for (int i = 0; i < hd.nframe; i++)
	{
		fmt::print(fp, "{:12.8f} ", hd.dt * i);
		for (auto&& ndx : allIndex)
		{
			if (contain(allPath[i].first, ndx)) // in center
			{
				fmt::print(fp, "{:5d} ", 1);
			}
			else if (contain(allPath[i].second, ndx)) // in surface
			{
				fmt::print(fp, "{:5d} ", 2);
			}
			else // in outer
			{
				fmt::print(fp, "{:5d} ", 0);
			}
		}
		for (int j = 0; j < allIndex.size(); j++)
		{
			fmt::print(fp, "{:5d} ", pairNumbr(i, j));
		}
		fmt::print(fp, "\n");
	}
	
	return 0;
}
