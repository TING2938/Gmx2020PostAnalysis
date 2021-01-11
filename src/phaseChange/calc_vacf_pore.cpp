/*  Author     : TING
 *  Date       : 2019/09/09
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculation vacf
 */
#include <itp/gmx>
 
struct Handle : public itp::GmxHandle
{
	using GmxHandle::GmxHandle;

	vecd calcVacf(const itp::Vector<matd>& vel, veci& index, int lateral)
	{
		double tmp;
		size_t halfFrame = vel.size() / 2;
		vecd vacf(halfFrame);
		vacf.fill(0);

		for (auto&&i : index)
		{
			for (size_t j = 0; j != halfFrame; ++j)
			{
				for (size_t k = 0; k != halfFrame; ++k)
				{
					tmp = 0;
					for (size_t m = 0; m != 3; ++m)
					{
						if (m != lateral)
						{
							tmp += vel[k + j][i][m] * vel[k][i][m];
						}
					}
					vacf[j] += tmp;
				}
			}
		}
		return vacf;
	}
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	hd.flags = TRX_READ_X | TRX_READ_V;
	hd.ngrps = 1;

	int lateral = 2;
	int dim1 = 0;
	real lowPos1 = -2; // nm;
	real upPos1 = 30; // nm;
	int dim2 = 2;
	real lowPos2 = -2;
	real upPos2 = 30;

	hd.pa = {
		{ "-lateral", FALSE, etINT, {&lateral}, "lateral to calc. 0(x), 1(y), 2(z)"},
		{ "-d1", FALSE, etINT, {&dim1}, "direction to calc. 0(x), 1(y), 2(z)"},
		{ "-up1", FALSE, etREAL, {&upPos1}, "up position of region of molecule/ion (nm)" },
		{ "-low1", FALSE, etREAL, {&lowPos1}, "low position of region of molecule/ion (nm)" },
		{ "-d2", FALSE, etINT, {&dim2}, "direction to calc. 0(x), 1(y), 2(z)"},
		{ "-up2", FALSE, etREAL, {&upPos2}, "up position of region of molecule/ion (nm)" },
		{ "-low2", FALSE, etREAL, {&lowPos2}, "low position of region of molecule/ion (nm)" }
	};

	hd.fnm = {
		{ efXVG, "-o", "vacf", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();

	auto velc = hd.initPosc(0);
	auto posc = hd.initPosc(0);

	itp::Vector<matd> vel;

	veci index0, index1;
	hd.loadPositionCenter(posc, 0);

	for (int i = 0; i != posc.nrow(); ++i)
	{
		if (lowPos1 < posc[i][dim1] && posc[i][dim1] < upPos1 && lowPos2 < posc[i][dim2] && posc[i][dim2] < upPos2)
		{
			index0.append(i);
		}
	}

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadVelocityCenter(velc, 0);
		hd.loadPositionCenter(posc, 0);
		vel.append(velc);

		index1.clear(); 
		for (int i = 0; i != posc.nrow(); ++i)
		{
			if (lowPos1 < posc[i][dim1] && posc[i][dim1] < upPos1 && lowPos2 < posc[i][dim2] && posc[i][dim2] < upPos2 && index0.contains(i))
			{
				index1.append(i);
			}
		} 
		index0 = index1;

	} while (hd.readNextFrame());

	auto vacf = hd.calcVacf(vel, index1, lateral);

	auto halfIndex = vel.size() / 2;
	for (int i = 0; i < halfIndex; i++)
	{
		vacf[i] = vacf[i] / halfIndex / index1.size() * 10000;
	}

	vecd normalvacf(halfIndex);
	for (int i = 0; i < halfIndex; i++)
	{
		normalvacf[i] = vacf[i] / vacf[0];
	}

	auto outputfileR = hd.openWrite(hd.get_ftp2fn(efXVG));
	for (int i = 0; i < halfIndex; i++)
	{
		fprintf(outputfileR, "%8.2d  %8.4f  %8.4f\n", i , vacf[i], normalvacf[i]);
	}
	fclose(outputfileR);

	return 0;
}

