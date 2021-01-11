#include <itp/gmx>

gmx_main(test)
{
	itp::GmxHandle hd(argc, argv); 
	hd.desc = {
		"calculate the temperature distrubution of group\n",
		"nbin: default 10."
	};
	hd.flags = TRX_READ_X | TRX_READ_V; 
	hd.ngrps = 1;  
	int nbin = 10;

	hd.pa = {
		{ "-nbin",  FALSE, etINT, {&nbin}, "nbins."}
	};

	hd.fnm = {
		{ efXVG, "-o", "Temperature", ffWRITE }
	};

	hd.init();

	hd.readFirstFrame();

	
	auto fp = itp::gmx_open(hd.get_opt2fn("-o"), argc, argv);
	/****************** do some analysis     *****************************************/

	vecd temp(nbin, 0);
	int frameNumber = 0;
	double boxZLength = 0.0;
	vecd Etotol(nbin, 0);
	vecu num(nbin, 0);
	auto nr = hd.ngx[0];

	do {

		boxZLength = hd.fr.box[2][2];
		Etotol.fill(0);
		num.fill(0);
		
		for (size_t n = 0; n != hd.ngrps; n++) {
			for (size_t i = 0; i != nr; i++) {
				auto ndx = hd.index[n][i];
				auto x = hd.fr.x[ndx];
				auto v = hd.fr.v[ndx];
				auto m = hd.top->atoms.atom[ndx].m;

				size_t index_z = static_cast<size_t>((x[2] <= boxZLength) ? (nbin * x[2] / boxZLength) : (nbin - 1));
				double Etmp = 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
				Etotol[index_z] += Etmp;
				++num[index_z];
			}
		}

		for (size_t i = 0; i != nbin; ++i) {
			double tmpT = Etotol[i] * 166.0539 / (1.5 * 1.3806505 * (num[i] - 1));
			temp[i] += tmpT;
		}

		++frameNumber;
	} while (hd.readNextFrame());

	double partLength = boxZLength / nbin;
	double avgTemperature = 0.0;
	for (size_t i = 0; i != nbin; ++i) {
		temp[i] /= frameNumber;
		avgTemperature += temp[i];
		fprintf(fp, "%10.3f", (2 * i + 1.0) * partLength / 2.0);
		fprintf(fp, "%10.3f\n", temp[i]);
	}

	fprintf(stderr, "The average Temperature of %s is %5.3f\n", hd.grpname[0], avgTemperature / nbin);
	/*********************************************************************************/
	fclose(fp);
	return 0;
}


