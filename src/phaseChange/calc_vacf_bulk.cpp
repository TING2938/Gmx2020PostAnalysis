/*  Author     : TING
 *  Date       : 2019/09/09
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculation vacf
 */
#include <itp/gmx>

 
class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	vecd calcVacf(const itp::Vector<matd>& vel, int grp)
	{
		double tmp;
		size_t halfFrame = vel.size() / 2;
		vecd vacf(halfFrame);
		vacf.fill(0);

		for (size_t i = 0; i != natoms[grp].size(); ++i)
		{
			for (size_t j = 0; j != halfFrame; ++j)
			{
				for (size_t k = 0; k != halfFrame; ++k)
				{
					tmp = 0;
					for (size_t m = 0; m != 3; ++m)
					{
						tmp += vel[k + j][i][m] * vel[k][i][m];
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

	hd.fnm = {
		{ efXVG, "-o", "vacf", ffWRITE }
	};

	hd.init();
	hd.readFirstFrame();

	auto velc = hd.initPosc(0);
	itp::Vector<matd> vel;

	/* ---------------------------------------------------------------------------- */
	do
	{
		hd.loadVelocityCenter(velc, 0);
		vel.append(velc);

	} while (hd.readNextFrame());

	auto vacf = hd.calcVacf(vel, 0);

	auto halfIndex = vel.size() / 2;
	for (int i = 0; i < halfIndex; i++)
	{
		vacf[i] = vacf[i] / halfIndex / hd.natoms[0].size() * 10000;
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

