/*  Author     : TING
 *  Date       : 2019/09/09
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculation ecacf
 */
#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;

	void calcCacf() 
	{
		size_t halfFrame = qv.size() / 2; 
		double tmp; 
		cacf.resize(halfFrame);
		cacf.fill(0); 

		for (size_t j = 0; j != halfFrame; ++j)
		{
			for (size_t k = 0; k != halfFrame; ++k)
			{
				tmp = 0;
				for (size_t m = 0; m != 3; ++m)
				{
					tmp += qv[k + j][m] * qv[k][m];
				}
				cacf[j] += tmp;
			}
		}
	} 

public:
	std::vector<std::array<double, 3>> qv;
	Eigen::ArrayXd cacf;
};

gmx_main(temp)
{
	double kb = 1.38064852e-23; // J/K

	Handle hd(argc, argv);

	hd.flags = TRX_READ_X | TRX_READ_V;
	hd.ngrps = 1;

	real temperature = 300; // K
	hd.pa = {
		{"-temp", FALSE, etREAL, {&temperature}, "system temperature (K)"}
	};

	hd.fnm = {
		{ efXVG, "-o", "IR", ffWRITE }
	};

	hd.init(false);
	hd.readFirstFrame();

	double qvx, qvy, qvz, q;
	double dt = 0, tmpt = 0;
	double volume = 0.0;
	
	int ndx;
	/* ---------------------------------------------------------------------------- */
	do
	{
		qvx = 0;
		qvy = 0;
		qvz = 0;
		for (int i = 0; i != hd.ngx[0]; ++i)
		{
			ndx = hd.index[0][i];
			q = hd.top->atoms.atom[ndx].q;
			qvx += q * hd.fr->v[ndx][XX];
			qvy += q * hd.fr->v[ndx][YY];
			qvz += q * hd.fr->v[ndx][ZZ];
		}
		hd.qv.push_back({ qvx, qvy, qvz });

		dt = tmpt;
		tmpt = hd.fr->time;
		dt = tmpt - dt;

		volume += (hd.Lbox[XX] * hd.Lbox[YY] * hd.Lbox[ZZ]);

	} while (hd.readNextFrame());

	volume /= hd.nframe;

	hd.calcCacf();

	auto halfIndex = hd.qv.size() / 2;
	Eigen::ArrayXd time(halfIndex);
	for (int i = 0; i != halfIndex; i++)
	{
		hd.cacf[i] /= (halfIndex);
		time[i] = i * dt;
	}
	
	Eigen::ArrayXd inteCacf = itp::cumtrapz(time, hd.cacf);
	Eigen::ArrayXd conductivity = inteCacf * (1 / (3 * volume * kb * temperature) * 2.56697e-17); // S/m

	auto file = hd.openWrite(hd.get_ftp2fn(efXVG));
	fmt::print(file, "# System volume is {:5.3f} nm^3\n", volume);
	fmt::print(file, "# System temperature is {:5.3f} K\n", temperature);
	fmt::print(file, "# time  cacf  conductivity (ps, e^2 nm^2 ps^-2, S/m)\n");

	for (int i = 0; i < halfIndex; i++)
	{
		fmt::print(file, "{:12.5f}  {:12.5f} {:12.5f}\n", time[i] , hd.cacf[i], conductivity[i]);
	}
	fclose(file);

	return 0;
}

