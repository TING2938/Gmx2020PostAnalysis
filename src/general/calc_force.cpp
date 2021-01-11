/*  Author     : TING
 *  Date       : 2019/05/05
 *  Email      : yeting2938@hust.edu.cn
 *  Desc       : calculate paired and unpair number along z axis in channel simulation.
 */
#include <itp/gmx>

class Handle : public itp::GmxHandle
{
public:
	using GmxHandle::GmxHandle;
};

gmx_main(temp)
{
	Handle hd(argc, argv);

	hd.ngrps = 1; /* number of group(s), grp1 and grp2, anion and cation. */
	hd.flags = TRX_READ_X | TRX_READ_F;
	/* add some user-defined pargs. */

	hd.fnm = {
		{ efXVG, "-o", "force", ffWRITE }
	};

	hd.init(false);
	hd.readFirstFrame();

	std::vector<double> force;
	double fx, fy, fz;

	/* ---------------------------------------------------------------------------- */
	do
	{
		fx = hd.fr->f[10800][0];
		fy = hd.fr->f[10800][1];
		fz = hd.fr->f[10800][2];
		force.push_back(std::sqrt(fx * fx + fy * fy + fz * fz));

	} while (hd.readNextFrame());

	auto file = hd.openWrite("force.xvg");
	for (auto&& i : force)
	{
		fmt::print(file, "{}\n", i);
	}
	fclose(file);

	return 0;
}
