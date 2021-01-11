/* yeting
 * 2020.12.14
 * yeting2938@hust.edu.cn
 * for onfly 3D density post ananysis
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <itp/getopt>
#include <iomanip>

using real = float;

struct OnflyPostHandle
{
public:
	OnflyPostHandle(std::string fnm);
	void readHead();
	bool readNextFrame();

public:
	std::ifstream fid;
	int headSize;
	std::string commandLine;
	int ngrps;
	int nRegion;
	int userint1;
	int userint2;
	double dt;
	std::vector<real> lowPos;
	std::vector<real> upPos;
	std::vector<real> dbin;
	std::vector<int> nbin;
	int totNbin;
	real Lbox[3];
	std::vector<std::vector<real>> massdens;
	std::vector<std::vector<real>> numdensM;
	std::vector<std::vector<real>> chargedens;
	std::vector<std::vector<real>> numdensG;
};

int main(int argc, char** argv)
{
	std::string inputFnm = "densMNC3D.onfly";
	std::string outputFnm = "onfly3D.dat";
	double beginTime = 0;
	double endTime = 0;

	itp::Getopt getopt(argc, argv);
	getopt.getFixPos(inputFnm, 1, true, "input file name");
	getopt(beginTime, "-b", false, "begin time (ps)");
	getopt(endTime, "-e", false, "end time (ps)");
	getopt(outputFnm, "-o", false, "output file name");
	getopt.finish();

	printf("beginTime: %f\n", beginTime);
	printf("endTine: %f\n", endTime);
	printf("input file name: %s\n", inputFnm.c_str());
	printf("output file name: %s\n", outputFnm.c_str());

	OnflyPostHandle hd(inputFnm);
	hd.readHead();

	printf("freq: %d\n", hd.userint2);
	printf("dt: %f\n", hd.dt);

	double dt = hd.dt * hd.userint2;
	size_t beginStep = beginTime / dt;
	size_t endStep = endTime / dt;
	if (endStep == 0)
	{
		endStep = -1;
	}

	std::vector<std::vector<real>> massdens;
	std::vector<std::vector<real>> numdensM;
	std::vector<std::vector<real>> chargedens;
	std::vector<std::vector<real>> numdensG;
	
	massdens.resize(hd.ngrps);
	chargedens.resize(hd.ngrps);
	numdensM.resize(hd.ngrps);
	numdensG.resize(hd.ngrps);
	for (int i = 0; i < hd.ngrps; i++)
	{
		massdens[i].resize(hd.totNbin, 0);
		numdensM[i].resize(hd.totNbin, 0);
		chargedens[i].resize(hd.totNbin, 0);
		numdensG[i].resize(hd.totNbin, 0);
	}

	for (size_t i = 0; i < beginStep; i++)
	{
		if (!hd.readNextFrame())
		{
			printf("Error, no more frame!\n");
			std::exit(-1);
		}
	}

	printf("Analysis ...\n");
	size_t totFrame = 0;
	for (size_t i = beginStep; i <= endStep; i++)
	{
		if (!hd.readNextFrame())
			break;
		for (int i = 0; i < hd.ngrps; i++)
		{
			for (int j = 0; j < hd.totNbin; j++)
			{
				massdens[i][j] += hd.massdens[i][j];
				numdensM[i][j] += hd.numdensM[i][j];
				chargedens[i][j] += hd.chargedens[i][j];
				numdensG[i][j] += hd.numdensG[i][j];
			}
		}
		totFrame++;
	}

	printf("Read frame from %.3f ps to %.3f ps, %zd frames are averaged\n", beginStep * dt, (beginStep + totFrame) * dt, totFrame);

	for (int i = 0; i < hd.ngrps; i++)
	{
		for (int j = 0; j < hd.totNbin; j++)
		{
			massdens[i][j] /= totFrame;
			numdensM[i][j] /= totFrame;
			chargedens[i][j] /= totFrame;
			numdensG[i][j] /= totFrame;
		}
	}

	std::ofstream ofile(outputFnm);

	ofile << "% onflyCommandline = \"" << hd.commandLine << "\"\n";
	ofile << "% Lbox = [";
	for (auto&& i : hd.Lbox)
		ofile << i << " ";
	ofile << "];\n";
	ofile << "% lowPos = [";
	for (auto&& i : hd.lowPos)
		ofile << i << " ";
	ofile << "];\n";
	ofile << "% upPos = [";
	for (auto&& i : hd.upPos)
		ofile << i << " ";
	ofile << "];\n";
	ofile << "% fnm = \"" << outputFnm << "\";\n";
	ofile << "% nbin = [";
	for (auto&& i : hd.nbin)
	{
		ofile << i << " ";
	}
	ofile << "];\n";
	ofile << "% nRegion = " << hd.nRegion << ";\n";
	ofile << "% [data, xyzRange, nbin, lowPos, upPos, Lbox] = loadOnflyData3D(fnm);\n";
	ofile << "%\n";
	ofile <<
		R"(% function [ret, xyzRange, nbin, lowPos, upPos, Lbox] = loadOnflyData3D(fnm)
% allData = importdata(fnm);
% data = allData.data;
% textData = allData.textdata;
% Lbox = double(regexp(string(textData{2}), "[\d.+-]+", "match"));
% lowPos = double(regexp(string(textData{3}), "[\d.+-]+", "match"));
% upPos = double(regexp(string(textData{4}), "[\d.+-]+", "match"));
% nbin = double(regexp(string(textData{6}), "[\d.+-]+", "match"));
% xyzRange = cell(1);
% for ii=1:length(nbin)
%    xyzRange{ii} = linspace(lowPos(ii), upPos(ii), nbin(ii)); 
% end
% [~, col] = size(data);
% nRegion = length(nbin) / 3;
% ret = cell(1);
% tmpIndex = 0;
% for ii = 1:nRegion
%     eachNbin = nbin(3*(ii-1)+1) * nbin(3*(ii-1)+2) * nbin(3*ii);
%     ret{ii} = squeeze(reshape(data((tmpIndex+1):(tmpIndex+eachNbin), :), [nbin((3*(ii-1)+1):3*ii), col]));
%     tmpIndex = tmpIndex + eachNbin;
% end
% if nRegion == 1
%     ret = ret{1};
% end
% end)";
	ofile << "\n";

	ofile.setf(std::ios::scientific);
	ofile.precision(8);
	for (int j = 0; j < hd.totNbin; j++)
	{
		for (int i = 0; i < hd.ngrps; i++)
		{
			ofile << std::setw(15) << massdens[i][j] << "\t" 
				<< std::setw(15) << numdensM[i][j] << "\t" 
				<< std::setw(15) << chargedens[i][j] << "\t" 
				<< std::setw(15) << numdensG[i][j] << "\t";
		}
		ofile << "\n";
	}
	ofile.close();
	printf("Analysis done.");
}

OnflyPostHandle::OnflyPostHandle(std::string fnm)
{
	fid.open(fnm, std::ios::binary);
	if (!fid)
	{
		printf("Cannot open file \"%s\", exit!\n", fnm.c_str());
		std::exit(-1);
	}
}

void OnflyPostHandle::readHead()
{
	int magicNumber = 0;
	fid.read((char*)&magicNumber, sizeof(int));
	if (magicNumber == 20080513)
	{
		fid.read((char*)&headSize, sizeof(int));
		int commandLineSize = 0;
		fid.read((char*)&commandLineSize, sizeof(int));
		commandLine.resize(commandLineSize);
		fid.read((char*)commandLine.data(), commandLine.size());
		fid.read((char*)&ngrps, sizeof(int));
		fid.read((char*)&nRegion, sizeof(int));
		fid.read((char*)&userint1, sizeof(int));
		fid.read((char*)&userint2, sizeof(int));
		fid.read((char*)&dt, sizeof(double));
		lowPos.resize(3 * nRegion);
		upPos.resize(3 * nRegion);
		nbin.resize(3 * nRegion);
		dbin.resize(3 * nRegion);
		fid.read((char*)lowPos.data(), sizeof(real) * 3 * nRegion);
		fid.read((char*)upPos.data(), sizeof(real) * 3 * nRegion);
		fid.read((char*)dbin.data(), sizeof(real) * 3 * nRegion);
		fid.read((char*)nbin.data(), sizeof(int) * 3 * nRegion);
		fid.read((char*)Lbox, sizeof(real) * 3);
	}
	else
	{
		std::cout << "Error format, cannot read head!" << std::endl;
		std::exit(-1);
	}
	totNbin = 0;
	for (int i = 0; i < nRegion; i++)
		totNbin += nbin[i * 3] * nbin[i * 3 + 1] * nbin[i * 3 + 2];
	massdens.resize(ngrps);
	chargedens.resize(ngrps);
	numdensM.resize(ngrps);
	numdensG.resize(ngrps);
	for (int i = 0; i < ngrps; i++)
	{
		massdens[i].resize(totNbin);
		chargedens[i].resize(totNbin);
		numdensM[i].resize(totNbin);
		numdensG[i].resize(totNbin);
	}
}

bool OnflyPostHandle::readNextFrame()
{
	int magicNumber = 0;
	while (true)
	{
		fid.read((char*)&magicNumber, sizeof(int));
		if (magicNumber == 20201210)
		{
			for (int i = 0; i < ngrps; i++)
			{
				fid.read((char*)(massdens[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(numdensM[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(chargedens[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(numdensG[i].data()), sizeof(real) * totNbin);
			}
			return true;
		}
		else if (magicNumber == 20080513)
		{
			int ignoreSize = 0;
			fid.read((char*)&ignoreSize, sizeof(int));
			fid.seekg(ignoreSize - 2 * sizeof(int), std::ios::cur);
			magicNumber = 0;
			continue;
		}
		else
		{
			return false;
		}
	}
}