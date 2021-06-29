#include <string>
#include <vector>
#include <map>
#include <array>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <itp/timer>


class NBCellSearch_Cutoff
{
    struct ElectrodeCell
    {
        int cellIndex; // index of this cell
        std::vector<int> atomIndex;   // index of atom in this cell
        std::vector<int> JcellIndex; // electrolyte cell index around this cell
    };

    struct ElectrolyteCell
    {
        std::vector<int> atomIndex;  // index of atom in this cell
    };

private:
    int Nelectrode;
    int NAtoms;
    std::vector<std::array<double, 3>> x;
    std::array<double, 3> csize, Lbox, halfLbox;

    double rcutoff, rcutoff2;
    std::array<int, 3> cellShape;
    std::vector<ElectrodeCell> electrodeCells;
    std::vector<ElectrolyteCell> electrolyteCells;
    std::vector<std::vector<int>> nblists, nblists_normal;

private:
    inline int index3to1(int i, int j, int k)
    {
        return i + j * cellShape[0] + k * cellShape[0] * cellShape[1];
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

    inline bool checkWithinCutoff(const std::array<int, 3>& a1, const std::array<int, 3>& a2)
    {
        double low1, low2, width1, width2, dcenter;
        double dx, sqDistance = 0;
        for (int m = 0; m < 3; m++)
        {
            low1 = a1[m] * csize[m];
            low2 = a2[m] * csize[m];
            width1 = (low1 + csize[m] > Lbox[m]) ? (Lbox[m] - low1) : csize[m];
            width2 = (low2 + csize[m] > Lbox[m]) ? (Lbox[m] - low2) : csize[m];
            dcenter = pbc_dx((low1 + width1 / 2.0) - (low2 + width2 / 2.0), m);
            dx = dcenter - (width1 + width2) / 2.0;
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
        electrolyteCells.resize(cellShape[0] * cellShape[1] * cellShape[2]);
    }

    void setCutoff(double cutoff)
    {
        rcutoff = cutoff;
        rcutoff2 = cutoff * cutoff;
    }

    void setNelectrode(int Nele)
    {
        Nelectrode = Nele;
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

    void calcElectrodeCells()
    {
        for (auto&& i : electrodeCells)
            i.atomIndex.clear();

        std::map<int, std::map<int, std::map<int, std::vector<int>>>> allMap;
        std::array<int, 3> me;
        double tmpx;

        for (int i = 0; i < Nelectrode; i++)
        {
            if (i == 1410)
            {
                int adf = 2;
            }
            for (int m = 0; m < 3; m++)
            {
                tmpx = x[i][m];
                if (tmpx < 0)
                    tmpx += Lbox[m];
                if (tmpx >= Lbox[m])
                    tmpx -= Lbox[m];
                me[m] = tmpx / csize[m];
            }
            allMap[me[0]][me[1]][me[2]].push_back(i);
        }

        for (auto&& m1 : allMap)
        {
            me[0] = m1.first;
            for (auto&& m2 : m1.second)
            {
                me[1] = m2.first;
                for (auto&& m3 : m2.second)
                {
                    me[2] = m3.first;
                    ElectrodeCell electrodeCell;
                    electrodeCell.cellIndex = index3to1(me);
                    m3.second.shrink_to_fit();
                    electrodeCell.atomIndex = std::move(m3.second);

                    for (int i = 0; i < cellShape[0]; i++)
                        for (int j = 0; j < cellShape[1]; j++)
                            for (int k = 0; k < cellShape[2]; k++)
                            {
                                if (checkWithinCutoff(me, { i, j, k }))
                                    electrodeCell.JcellIndex.push_back(index3to1(i, j, k));
                            }
                    std::sort(electrodeCell.JcellIndex.begin(), electrodeCell.JcellIndex.end());
                    electrodeCell.JcellIndex.shrink_to_fit();
                    electrodeCells.push_back(std::move(electrodeCell));
                }
            }
        }
        electrodeCells.shrink_to_fit();
    }

    void calcElectrolyteCells()
    {
        for (auto&& i : electrolyteCells)
            i.atomIndex.clear();

        double tmpx;
        std::array<int, 3> me;
        for (int i = Nelectrode; i < NAtoms; i++)
        {
            if (i == 12604)
            {
                int adfd = 1;
            }
            for (int m = 0; m < 3; m++)
            {
                tmpx = x[i][m];
                if (tmpx < 0)
                    tmpx += Lbox[m];
                if (tmpx >= Lbox[m])
                    tmpx -= Lbox[m];
                me[m] = tmpx / csize[m];
            }
            electrolyteCells[index3to1(me)].atomIndex.push_back(i);
        }
    }

    void calcNBlists()
    {
        for (auto&& nbl : nblists)
        {
            nbl.clear();
        }
        std::array<double, 3> dx;
        double r2;

        for (int i = 0; i < electrodeCells.size(); i++)
        {
            for (auto&& p : electrodeCells[i].atomIndex)
            {
                for (auto&& j : electrodeCells[i].JcellIndex)
                {
                    for (auto&& q : electrolyteCells[j].atomIndex)
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
                                nblists[p].push_back(q);
                            }
                        }
                    }
                }
            }
        }
        for (auto&& nbl : nblists)
        {
            std::sort(nbl.begin(), nbl.end());
        }
    }

    void calcNBlists_normal()
    {
        for (auto&& nbl : nblists_normal)
            nbl.clear();

        std::array<double, 3> dx;
        double r2;
        double rcutoff2 = rcutoff * rcutoff;

        for (int p = 0; p < Nelectrode; p++)
        {
            for (int q = Nelectrode; q < NAtoms; q++)
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
                        nblists_normal[p].push_back(q);
                    }
                }
            }
        }
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
        for (int i = 0; i < 100; i++)
        {
            double cutoff = (0.2 + i * 0.02) * 1.2;
            NBCellSearch_Cutoff nbs;
            nbs.readGro();
            nbs.setCellSize({ cutoff*1.2, cutoff, cutoff*1.2 });
            nbs.setNelectrode(10800);
            nbs.setCutoff(cutoff);
            nbs.calcElectrodeCells();

            nbs.calcElectrolyteCells();
            nbs.calcNBlists();

            nbs.calcNBlists_normal();

            if (!nbs.compare())
            {
                printf("not equal!, i=%d\n", i);
                std::exit(-1);
            }
        }
    }
    else
    {
        itp::Timeit timeit[3];

        double cutoff = 1.2 * 1.2;
        NBCellSearch_Cutoff nbs;
        nbs.readGro();
        nbs.setCellSize({ cutoff * 0.4, cutoff * 0.4, cutoff * 0.4 });
        nbs.setNelectrode(10800);
        nbs.setCutoff(cutoff);
        nbs.calcElectrodeCells();

        for (int i = 0; i < 10; i++)
        {
            timeit[0].start();
            nbs.calcElectrolyteCells();
            timeit[0].pause();

            timeit[1].start();
            nbs.calcNBlists();
            timeit[1].pause();

            timeit[2].start();
            nbs.calcNBlists_normal();
            timeit[2].pause();
        }

        timeit[0].printSpan("cell elety: ", "s\n");
        timeit[1].printSpan("cell nblist: ", "s\n");
        timeit[2].printSpan("normol: ", "s\n");
        if (!nbs.compare())
        {
            printf("not equal!\n");
        }
    }
}