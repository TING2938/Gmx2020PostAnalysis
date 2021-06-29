#include <string>
#include <vector>
#include <map>
#include <array>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <itp/timer>
#include <omp.h>

class NBCellSearch_Cutoff
{
private:
    int Nelectrode;
    int NAtoms;
    std::vector<std::array<double, 3>> x;
    std::array<double, 3> csize, Lbox, halfLbox;

    double rcutoff, rcutoff2;
    std::array<int, 3> cellShape;
    int cellShapeXY;
    std::vector<std::vector<int>> electrodeCellLists;  // cell index of electrolyte cell around each electrode atom
    std::vector<std::vector<int>> electrolyteCellAtomLists; // atom index lists in each electrolyte cell 
    std::vector<std::vector<int>> nblists, nblists_normal;

private:
    inline int index3to1(int i, int j, int k)
    {
        return i + j * cellShape[0] + k * cellShapeXY;
    }
    inline int index3to1(const std::array<int, 3>& index3)
    {
        return index3to1(index3[0], index3[1], index3[2]);
    }

    inline double pbc_dx(double dx, int dim)
    {
        dx = std::abs(dx);
        if (dx > Lbox[dim] / 2.0)
            return Lbox[dim] - dx;
        return dx;
    }

    inline bool checkWithinCutoff(const std::array<double, 3>& coordinate, const std::array<int, 3>& cellIndex)
    {
        double low, width;
        double dx, sqDistance = 0;
        for (int m = 0; m < 3; m++)
        {
            low = cellIndex[m] * csize[m];
            width = (low + csize[m] > Lbox[m]) ? (Lbox[m] - low) : csize[m];
            dx = pbc_dx((low + width / 2.0) - coordinate[m], m) - width / 2.0;
            if (dx > 0)
            {
                sqDistance += dx * dx;
                if (sqDistance >= rcutoff2)
                    return false;
            }
        }
        return true;
    }

public:
    void setCellSize(const std::array<double, 3>& size)
    {
        csize = size;
        for (int i = 0; i < 3; i++)
        {
            cellShape[i] = int(std::ceil(Lbox[i] / csize[i]));
        }
        electrolyteCellAtomLists.resize(cellShape[0] * cellShape[1] * cellShape[2]);
        cellShapeXY = cellShape[0] * cellShape[1];
    }

    void setCutoff(double cutoff)
    {
        rcutoff = cutoff;
        rcutoff2 = cutoff * cutoff;
    }

    void setNelectrode(int Nele)
    {
        Nelectrode = Nele;
        electrodeCellLists.resize(Nelectrode);
        nblists.resize(Nelectrode);
        nblists_normal.resize(Nelectrode);
        for (int i = 0; i < Nelectrode; i++)
        {
            nblists[i].reserve((NAtoms - Nelectrode) / 2);
            nblists_normal[i].reserve((NAtoms - Nelectrode) / 2);
        }
    }

    void readGro()
    {
        // load .gro file;
        std::ifstream file("mofs.gro");
        if (!file)
        {
            printf("Cannot open file \"mofs.gro\"!\n");
            std::exit(-1);
        }
        std::stringstream ss;
        std::string line;

        std::getline(file, line);
        std::getline(file, line);

        ss.str(line);
        ss >> NAtoms;
        ss.clear();

        x.resize(NAtoms);
        for (int i = 0; i < NAtoms; i++) {
            std::getline(file, line);
            ss.str(line.substr(20, 24));
            ss >> x[i][0] >> x[i][1] >> x[i][2];
            ss.clear();
        }
        std::getline(file, line);
        ss.str(line);
        ss >> Lbox[0] >> Lbox[1] >> Lbox[2];
        for (int m = 0; m < 3; m++)
            halfLbox[m] = Lbox[m] / 2.0;
        ss.clear();
        file.close();
    }

    void calcElectrodeCellLists()
    {
        for (auto&& cell : electrodeCellLists)
            cell.clear();

        for (int i = 0; i < Nelectrode; i++)
        {
            for (int p = 0; p < cellShape[0]; p++)
                for (int q = 0; q < cellShape[1]; q++)
                    for (int r = 0; r < cellShape[2]; r++)
                    {
                        if (checkWithinCutoff(x[i], { p, q, r }))
                            electrodeCellLists[i].push_back(index3to1(p, q, r));
                    }
            electrodeCellLists[i].shrink_to_fit();
        }
    }

    void calcElectrolyteCellAtomLists()
    {
        for (auto&& atomList : electrolyteCellAtomLists)
            atomList.clear();

        double tmpx;
        std::array<int, 3> me;
        for (int i = Nelectrode; i < NAtoms; i++)
        {
            for (int m = 0; m < 3; m++)
            {
                tmpx = x[i][m];
                if (tmpx < 0)
                    tmpx += Lbox[m];
                if (tmpx >= Lbox[m])
                    tmpx -= Lbox[m];
                me[m] = tmpx / csize[m];
            }
            electrolyteCellAtomLists[index3to1(me)].push_back(i);
        }
    }

    void calcNBlists()
    {
        for (auto&& nbl : nblists)
        {
            nbl.clear();
        }
        
#pragma omp parallel for
        for (int i = 0; i < Nelectrode; i++)
        {
            for (auto&& cellIndex : electrodeCellLists[i])
            {
                for (auto&& j : electrolyteCellAtomLists[cellIndex])
                {
                    std::array<double, 3> dx;
                    dx[0] = pbc_dx(x[i][0] - x[j][0], 0);
                    if (dx[0] >= rcutoff)
                        continue;
                    dx[1] = pbc_dx(x[i][1] - x[j][1], 1);
                    if (dx[1] >= rcutoff)
                        continue;
                    dx[2] = pbc_dx(x[i][2] - x[j][2], 2);
                    if (dx[2] >= rcutoff)
                        continue;
                    double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                    if (r2 < rcutoff2)
                    {
                        nblists[i].push_back(j);
                    }
                }
            }
            // std::sort(nblists[i].begin(), nblists[i].end());
        }
    }

    void calcNBlists_normal()
    {
        for (auto&& nbl : nblists_normal)
            nbl.clear();

        double rcutoff2 = rcutoff * rcutoff;

#pragma omp parallel for 
        for (int p = 0; p < Nelectrode; p++)
        {
            for (int q = Nelectrode; q < NAtoms; q++)
            {
                std::array<double, 3> dx;
                for (int m = 0; m < 3; m++)
                {
                    dx[m] = pbc_dx(x[p][m] - x[q][m], m);
                }
                if (dx[0] < rcutoff && dx[1] < rcutoff && dx[2] < rcutoff)
                {
                    double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                    if (r2 < rcutoff2)
                    {
                        nblists_normal[p].push_back(q);
                    }
                }
            }
        }
    }

    std::vector<std::vector<int>> calcNBlists_normal_reverse()
    {
        std::vector<std::vector<int>> ret(NAtoms - Nelectrode);

        std::array<double, 3> dx;
        double r2;
        double rcutoff2 = rcutoff * rcutoff;

        for (int p = Nelectrode; p < NAtoms; p++)
        {
            for (int q = 0; q < Nelectrode; q++)
            {
                for (int m = 0; m < 3; m++)
                {
                    dx[m] = pbc_dx(x[p][m] - x[q][m], m);
                }
                if (dx[0] < rcutoff && dx[1] < rcutoff && dx[2] < rcutoff)
                {
                    r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                    if (r2 < rcutoff2)
                    {
                        ret[p-Nelectrode].push_back(q);
                    }
                }
            }
        }
        return ret;
    }

    std::vector<std::vector<int>> reverseNBlists()
    {
        std::vector<std::vector<int>> ret(NAtoms - Nelectrode);
        for (int i = 0; i < Nelectrode; i++)
        {
            for (auto&& j : nblists_normal[i])
            {
                ret[j - Nelectrode].push_back(i);
            }
        }
        return ret;
    }

    bool compare()
    {
        for (int i = 0; i < Nelectrode; i++)
        {
            if (nblists[i] != nblists_normal[i])
            {
                printf("in electrode: %d\n", i);
                return false;
            }
        }
        return true;
    }
};

int main()
{
    if (0)
    {
        for (int i = 0; i < 400; i++)
        {
            double cutoff = (0.2 + i * 0.005) * 1.2;
            NBCellSearch_Cutoff nbs;
            nbs.readGro();
            nbs.setCellSize({ cutoff*1.2, cutoff, cutoff*1.2 });
            nbs.setNelectrode(10800);
            nbs.setCutoff(cutoff);
            nbs.calcElectrodeCellLists();

            nbs.calcElectrolyteCellAtomLists();
            nbs.calcNBlists();

            nbs.calcNBlists_normal();

            if (!nbs.compare())
            {
                printf("not equal!, i=%d\n", i);
                std::exit(-1);
            }

            /*auto l1 = nbs.calcNBlists_normal_reverse();
            auto l2 = nbs.reverseNBlists();

            for (int i = 0; i < l1.size(); i++)
            {
                if (l1[i] != l2[i])
                {
                    printf("error of reverse, i=%d\n", i);
                    std::exit(-2);
                }
            }*/
        }
    }
    else
    {
        itp::Timeit timeit[3];

        std::vector<double> allTime(100-5+1);

        for (int i = 5; i <= 100; i++)
        {
            double cutoff = 1.2 * 1.2;
            NBCellSearch_Cutoff nbs;
            nbs.readGro();
            double rate = cutoff * i * 0.02;
            nbs.setCellSize({ rate, rate, rate });
            nbs.setNelectrode(10800);
            nbs.setCutoff(cutoff);
            nbs.calcElectrodeCellLists();

            timeit[0].start();
            nbs.calcElectrolyteCellAtomLists();
            nbs.calcNBlists();
            timeit[0].stop();

            allTime[i - 5] = timeit[0].span();
        }

        auto maxptr = std::max_element(allTime.begin(), allTime.end());
        auto minptr = std::min_element(allTime.begin(), allTime.end());

        printf("max[%d]: %f, min[%d]: %f\n", maxptr - allTime.begin(), *maxptr, minptr - allTime.begin(), *minptr);
    }
}

