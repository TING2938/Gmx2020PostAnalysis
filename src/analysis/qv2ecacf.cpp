#include <itp/core>
#include <itp/fileio>
#include <itp/utility>
#include <itp/getopt>
#include <fstream>
#include <omp.h>

int main(int argc, char** argv)
{
    double kb = 1.38064852e-23; // J/K
    double volume = 89; // nm^3
    double temperature = 298; // K
    std::string fnm = "qv.xvg";
    std::string outputFnm = "IR.xvg";
    int nbin = 100;
    double dt = 0.002; // ps
    double lowPos = 0;
    double upPos = 10;
    int nRegion = 2;
    std::vector<std::vector<double>> region;
    
    // region[0].pu ... 

    itp::Getopt getopt(argc, argv, "calc qv to ecacf");
    getopt(volume, "-v", true, "volume");
    getopt(temperature, "-t", true, "temperature");
    getopt(fnm, "-if", true, "input file name");
    getopt(outputFnm, "-of", true, "output file name");
    getopt(nbin, "-nbin", true, "nbin");
    getopt(dt, "-dt", true, "delta t");
    getopt(lowPos, "-low", true, "lowPos");
    getopt(upPos, "-up", true, "upPos");
    getopt(nRegion, "-nRegion", true, "nRegion");
    std::vector<double> tmpRegion;
    getopt.getArray(tmpRegion, "-region", 2*nRegion, true, "region");
    getopt.finish();
    region.resize(nRegion);
    for (int i = 0; i < nRegion; i++)
    {
        region[i].resize(2);
        region[i][0] = tmpRegion[2 * i];
        region[i][1] = tmpRegion[2 * i + 1];
    }

    std::vector<std::vector<int>> index;
    index.resize(nRegion);
    for (int i = 0; i < nRegion; i++)
    {
        index[i].resize(2);
        for (int m = 0; m < 2; m++)
        {
            index[i][m] = (region[i][m] - lowPos) / (upPos - lowPos) * nbin;
        }
    }
    int nthreads = omp_get_num_threads();

    auto data = itp::loadtxt(fnm);
    int nrow = data.rows();
    int nframe = nrow / nbin; // length of time

    std::vector<std::vector<Eigen::Array3d>> qv;
    qv.resize(nRegion);
    for (auto&& qvi : qv)
    {
        qvi.resize(nframe);
    }
    for (int n = 0; n < nRegion; n++)
    {
        int index_low = index[n][0];
        int index_up = index[n][1];

        #pragma omp parallel for num_threads(nthreads)
        for (int i = 0; i < nframe; i++)
        {
            qv[n][i] = data.middleRows(i * nbin + index_low, index_up - index_low).colwise().sum();
        }
    }

    int halframe = nframe / 2;
    Eigen::ArrayXXd cacf(nRegion, halframe);
    cacf.fill(0);


    for (int i = 0; i < nRegion; i++)
    {
    #pragma omp parallel for num_threads(nthreads)
        for (int j = 0; j < halframe; ++j)
        {
            for (int k = 0; k < halframe; ++k)
            {
                double tmp = 0;
                for (int m = 0; m < 3; ++m)
                {
                    tmp += qv[i][k + j][m] * qv[i][k][m];
                }
                cacf(i, j) += tmp;
            }
        }
    }

    Eigen::ArrayXd time(halframe);
    for (int i = 0; i < halframe; i++)
    {
        time[i] = i * dt;
    }
    for (int i = 0; i < nRegion; i++)
    {
        cacf.row(i) /= halframe;
    }

    Eigen::ArrayXXd inteCacf(nRegion, halframe);
    Eigen::ArrayXXd conductivity(nRegion, halframe);
    inteCacf.fill(0);
    conductivity.fill(0);
    for (int i = 0; i < nRegion; i++)
    {
        Eigen::ArrayXd row = cacf.row(i);
        inteCacf.row(i) = itp::cumtrapz(time, row);
        conductivity.row(i) = inteCacf.row(i) * (1 / (3 * volume * kb * temperature) * 2.56697e-17); // S/m
    }
    
    std::ofstream ofile(outputFnm);
    for (int i = 0; i < halframe; i++)
    {
        fmt::print(ofile, "{:12.5} ", time[i]);
        for (int n = 0; n < nRegion; n++)
        {
            fmt::print(ofile, "{:12.5} {:12.5}", cacf(n, i), conductivity(n, i));
        }
        fmt::print(ofile, "\n");
    }
    ofile.close();
}

