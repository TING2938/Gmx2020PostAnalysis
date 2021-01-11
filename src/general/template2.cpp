#include <cmath>
#include <gromacs/commandline/pargs.h>
#include <gromacs/fileio/trxio.h> 
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include <gromacs/topology/index.h>
#include <gromacs/topology/topology.h>
#include <gromacs/utility/arraysize.h>
#include <gromacs/utility/smalloc.h>
#include <gromacs/utility/fatalerror.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <gromacs/commandline/cmdlineinit.h>
#include <itp/core>

double* CreateVector(int rows)
{
	return new double[rows]();
}

int my_main(int argc, char* argv[])
{
	const char* desc[] = { 
		"calculate density distribution, only used for study."
	};

	static int   nbin = 300;
	static float lowPos = 0;
	static float upPos = 12;

	t_pargs pa[] = {
	  { "-nbin",FALSE, etINT, {&nbin},
		"number of bin sets for number"
	  },
	  { "-low",FALSE, etREAL, {&lowPos},
		"low position of region of molecule/ion (nm)"
	  },
	  { "-up",FALSE, etREAL, {&upPos},
		"up position of region of molecule/ion (nm)"
	  }
	};

	t_filenm fnm[] = {
		{ efTRX, "-f", NULL,  ffREAD },
		{ efNDX, NULL, NULL,  ffOPTRD },
		{ efTPR, nullptr, nullptr,  ffREAD },
		{ efXVG, "-o", "density", ffWRITE },
	};

	t_topology       *top;
	int              ePBC; 
	gmx_output_env_t *oenv;
	t_trxframe       fr;
	t_trxstatus      *status;
	int              flags = TRX_READ_X;   /* read only position */

	char             **grpname;
	static int       ngrps = 1;
	int              **index;
	int              *ngx;

#define NFILE asize(fnm)
	if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME,
		NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv)) {
		return 0;
	}

	top = read_top(ftp2fn(efTPR, NFILE, fnm), &ePBC); /* read topology file */
	snew(grpname, ngrps);
	snew(index, ngrps);
	snew(ngx, ngrps);

	auto aa = top->atoms.atom[0].resind;
	auto bb = top->atoms.resinfo[aa].nr;

	get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), ngrps, ngx, index, grpname);

	read_first_frame(oenv, &status, ftp2fn_null(efTRX, NFILE, fnm), &fr, flags);

	


    return 0;
}


int main(int argc, char *argv[])
{
    return gmx_run_cmain(argc, argv, &my_main);
}
