
#ifdef HAVE_CONFIG_H
#include "gmxpre-config.h"
#endif
#include <cmath>
#include <gromacs/commandline/pargs.h>
#include <gromacs/fileio/trxio.h>
#include <gromacs/fileio/xvgr.h>
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include <gromacs/topology/index.h>
#include <gromacs/topology/topology.h>
#include <gromacs/utility/arraysize.h>
#include <gromacs/utility/smalloc.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <gromacs/commandline/cmdlineinit.h>

int load_alpha_epsilon(int atom1, int atom2, double **A, double **E, double *a1, double *a2, double *e1, double *e2, int LJtype);
int load_positionCenter(t_trxframe fr, int nm, int atoms, int *index, double **pos, double Lbox[3], double *deloc, double totNCM);
int load_position(t_trxframe fr, int nm, int atoms, int *index, double ***posmolecule, double Lbox[3]);
double periodicity(double dx, double box);
double * CreateVector(int cols);
double **CreateMatrix(int rows, int cols);
double *** CreateMatrix_3d(int m, int n, int p)
{
    int i;
    double ***m2;
    m2 = new double**[m];
    for (i = 0; i < m; i++)
        m2[i] = CreateMatrix(n, p);
    return m2;
}


int my_main(int argc, char *argv[])
{
    const char *desc[] = {
      "this is a small test program meant to serve as a template ",
      "when writing your own analysis tools. The advantage of ",
      "using gromacs for this is that you have access to all ",
      "information in the topology, and your program will be ",
      "able to handle all types of coordinates and trajectory ",
      "files supported by gromacs. Go ahead and try it! ",
      "This test version just writes the coordinates of an ",
      "arbitrary atom to standard out for each frame. You can ",
      "select which atom you want to examine with the -n argument."
    };
    static int preframe = 0;
    static int endframe = 9999999;
    static int nStart1 = 0;        // starting index of reference atoms
    static int nStart2 = 10920;
    static int nm1 = 840;
    static int nm2 = 840;
    static int tN1 = 2;
    static int tN2 = 4;
    static int n_NCM = 2;
    static int DIMN = 2;
    static int LJtype = 1;
    static int   nbin = 300;
    static float D = 1.0;
    static float Cr = 0.0;
    static float CR = 0.8;
    static float lowPos = 0;
    static float upPos = 12;


    /* GMX-467 support data structure */
    static const char *normtype[] = { NULL, "no", "x", "y", "z", NULL };
    static const char *axtitle[] = { NULL, "no", "x", "y", "z", NULL };
    static gmx_bool    bTen = FALSE;
    static gmx_bool    bMW = TRUE;
    static gmx_bool    bRmCOMM = FALSE;
	
    t_pargs pa[] = {
      { "-tN1",FALSE, etINT, {&tN1},
        "the type of reference atom or molecule"
      },
      { "-tN2",FALSE, etINT, {&tN2},
        "the type of surrounding atom or molecule"
      },
      { "-NCM",FALSE, etINT, {&n_NCM},
        "center of molecule\n0: number; 1: charge; 2: mass"
      },
      { "-dim",FALSE, etINT, {&DIMN},
        "dim of bin sets\ndim>2 means getting D in 3-dimension,0 X,1 Y,2 Z"
      },
      { "-D",FALSE, etREAL, {&D},
        "dielectric constants"
      },
      { "-pre",FALSE, etINT, {&preframe},
        "Number of preframe"
      },
      { "-end",FALSE, etINT, {&endframe},
        "End number of trajectory set"
      },
      { "-Cr",FALSE, etREAL, {&Cr},
        "Minimum radius of spacing(nm)"
      },
      { "-CR",FALSE, etREAL, {&CR},
        "Maximam radius of spacing(nm)"
      },
      { "-nbin",FALSE, etINT, {&nbin},
        "number of bin sets for number"
      },
      { "-low",FALSE, etREAL, {&lowPos},
        "low position of region of molecule/ion (nm)"
      },
      { "-up",FALSE, etREAL, {&upPos},
        "up position of region of molecule/ion (nm)"
      },
      { "-type",FALSE, etINT, {&LJtype},
        "L-J type. 1:Lorentz-Berthelot; 2:geometric mean."
      }
    };

    /* GMX-467 support data structure */
    t_filenm           fnm[] = {
        { efTRX, "-f", NULL,  ffREAD },
        { efNDX, NULL, NULL,  ffOPTRD },
        { efTPR, nullptr, nullptr,  ffREAD },
        { efXVG, NULL, "msd", ffWRITE },
        { efXVG, "-mol", "diff_mol", ffOPTWR },
        { efPDB, "-pdb", "diff_mol", ffOPTWR }
    };




    /* GMX-467 support data structure */
    t_topology      *top;
    int             ePBC;
    matrix          box;
    char            title[256];
    const char      *trx_file, *tps_file, *ndx_file, *msd_file, *mol_file, *pdb_file;
    rvec            *xdum;
    gmx_bool        bTop;
    int             axis, type;
    real            dim_factor;
    gmx_output_env_t    *oenv;
    t_trxframe      fr;
    t_trxstatus     *status;
    int             flags = TRX_READ_X;   /* read only position */

    char            **grpname;
    static int      ngrps = 2;
    int         **indexxx;
    int			    *ngx;

    /* variables for this analysis */
    int     i, j;
    FILE    *inputPara, *outputfileR;
    int    NUM[999], numbin[999], T[99999], Number[999], numAt[10], meZ, index, g, l, h, step, k, atoms1, atoms2, atoms_r, nType, m, n, sameMol, nij = 0;
    double  atNCM[5][10][100], mol_NCM[3][10], totV[99999], totN[99999], totC[99999], Vv2[99999], Cc2[99999], dR[3], dr[3], Lbox[3], L[3], r, R, dRc, drc, dt, tmpt, dc, tV, tC, tN, tVC, dL[3];
    char    name0[60], name1[10], temp[36];
    double  ***pos1, ***pos2, **posc1, **posc2, **Atable, **Etable, **V, **C, **N, **Vv0, **Cc0;
    float	VDW[50][50], COU[50][50], Vv1[999], Cc1[999], dbin;//Vv1[900],tN1 number<900



#define NFILE asize(fnm)
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
        return 0;

    top = read_top(ftp2fn(efTPR, NFILE, fnm), &ePBC); /* read topology file */
    snew(grpname, ngrps);
    snew(indexxx, ngrps);
    snew(ngx, ngrps);

    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm), ngrps, ngx, indexxx, grpname);

    // get the atom information of geometry center
    // MSD_charge_para.dat comes from RDF_charge_para.dat
    if ((inputPara = fopen("MSD_charge_para_Energy.dat", "r")) == NULL) {
        printf("\nNOTE: Parameter file cannot be loaded, mission terminated.\n");
        exit(0);
    }
    fscanf(inputPara, "%s %d \n", temp, &nType);
    if (tN1 > nType) {
        printf("\nNOTE: there is NO type %d in MSD_charge_para_Energy.dat.\n", tN1);
        exit(0);
    }
    for (j = 0; j < nType; j++) {
        fscanf(inputPara, "%s %d \n", temp, &numAt[j]);
        for (i = 0; i < numAt[j]; i++) {
            fscanf(inputPara, "%d%s%d%s%s%lf%lf%lf%lf%lf ", &atoms_r, temp, &atoms_r, temp, temp, &atNCM[0][j][i], &atNCM[1][j][i], &atNCM[2][j][i], &atNCM[3][j][i], &atNCM[4][j][i]);
            if (atNCM[0][j][i] > 0) atNCM[0][j][i] = 1.0;
            else atNCM[0][j][i] = 0.0;
            //printf("         %d %d %.1f %8.4f %8.4f %8.4f %8.4f \n",j,nType,atNCM[0][j][i],atNCM[1][j][i],atNCM[2][j][i],atNCM[3][j][i],atNCM[4][j][i]);

        }
    }
    tN1--; tN2--;  // make it as index in array

    atoms1 = numAt[tN1]; atoms2 = numAt[tN2];
    // total mass/effective number/charge in a whole molecule (for solvent, the net charge is zero) 8/3/2010
    for (j = 0; j < nType; j++) {
        for (k = 0; k < 3; k++) {
            mol_NCM[k][j] = 0.0;
            for (i = 0; i < numAt[j]; i++)
                mol_NCM[k][j] += atNCM[k][j][i];

            if (abs(mol_NCM[k][j]) < 1e-6) mol_NCM[k][j] = 1.0;
        }
    }
    printf("\n# of atoms in reference type: %10d\n", atoms1);
    printf("Total:     %.1f %8.4f %8.4f \n", mol_NCM[0][tN1], mol_NCM[1][tN1], mol_NCM[2][tN1]);
    printf("           ---------------------\n");
    for (i = 0; i < atoms1; i++)
        printf("           %.1f %8.4f %8.4f %8.4f %8.4f \n", atNCM[0][tN1][i], atNCM[1][tN1][i], atNCM[2][tN1][i], atNCM[3][tN1][i], atNCM[4][tN1][i]);

    printf("\n# of atoms in surrounding type: %10d\n", atoms2);
    printf("Total:     %.1f %8.4f %8.4f \n", mol_NCM[0][tN2], mol_NCM[1][tN2], mol_NCM[2][tN2]);
    printf("           ---------------------\n");
    for (i = 0; i < atoms2; i++)
        printf("           %.1f %8.4f %8.4f %8.4f %8.4f \n", atNCM[0][tN2][i], atNCM[1][tN2][i], atNCM[2][tN2][i], atNCM[3][tN2][i], atNCM[4][tN2][i]);

    /******************************/

    nm1 = ngx[0] / atoms1;
    nm2 = ngx[1] / atoms2;
    printf("\n nm1 = %10d, nm2 = %10d\n", nm1, nm2);

    /******************************/

    /* GMX-467 support build-in function:
     * read_first_frame(oenv,&status,trx_file,&fr,flags);
     * read_next_frame(oenv,status,&fr);
     */
    read_first_frame(oenv, &status, ftp2fn_null(efTRX, NFILE, fnm), &fr, flags);
    for (step = 0; step < preframe; step++)
        read_next_frame(oenv, status, &fr);
    tmpt = fr.time;
    /* =================== Main body of code ================== */

    pos1 = CreateMatrix_3d(nm1, atoms1, 3);
    pos2 = CreateMatrix_3d(nm2, atoms2, 3);  // pos of atoms
    posc1 = CreateMatrix(nm1, 3);   // center of mols
    posc2 = CreateMatrix(nm2, 3);
    Atable = CreateMatrix(atoms1, atoms2);
    Etable = CreateMatrix(atoms1, atoms2);
    V = CreateMatrix(99999, nbin);
    C = CreateMatrix(99999, nbin);
    N = CreateMatrix(99999, nbin);
    Vv0 = CreateMatrix(nm1, nm2);
    Cc0 = CreateMatrix(nm1, nm2);
    index = 0;

    meZ = 0;

    do {
        Lbox[0] = fr.box[XX][XX]; Lbox[1] = fr.box[YY][YY]; Lbox[2] = fr.box[ZZ][ZZ];
        load_alpha_epsilon(atoms1, atoms2, Atable, Etable, atNCM[3][tN1], atNCM[3][tN2], atNCM[4][tN1], atNCM[4][tN2], LJtype);
        load_position(fr, nm1, atoms1, indexxx[0], pos1, Lbox);
        load_positionCenter(fr, nm1, atoms1, indexxx[0], posc1, Lbox, atNCM[n_NCM][tN1], mol_NCM[n_NCM][tN1]);
        /*if (tN1 == tN2) {
            for (i = 0; i < nm1; i++)
                for (j = 0; j < atoms1; j++)
                    for (k = 0; k < 3; k++)
                        pos2[i][j][k] = pos1[i][j][k];
            posc2[i][k] = posc1[i][k];
        }
        else*/
        load_position(fr, nm2, atoms2, indexxx[1], pos2, Lbox);
        load_positionCenter(fr, nm2, atoms2, indexxx[1], posc2, Lbox, atNCM[n_NCM][tN2], mol_NCM[n_NCM][tN2]);

        for (m = 0; m < atoms1; m++) { for (n = 0; n < atoms2; n++) { VDW[m][n] = 0.0; COU[m][n] = 0.0; } }
        for (i = 0; i < nbin; i++) { Vv2[i] = 0.0; Cc2[i] = 0.0; NUM[i] = 0.0; numbin[i] = 0; }

        dbin = 0; dbin = (upPos - lowPos) / nbin;
        for (i = 0; i < nm1; i++) {
            Number[i] = 0.0;	Vv1[i] = 0.0; Cc1[i] = 0.0;
            dL[DIMN] = 0;
            dL[DIMN] = posc1[i][DIMN];
            if (lowPos < dL[DIMN] && dL[DIMN] < upPos) {
                meZ = (dL[DIMN] - lowPos) / dbin; 
                numbin[meZ]++;
                for (j = 0; j < nm2; j++) {
                    dRc = 0.0; R = 0; Vv0[i][j] = 0.0; Cc0[i][j] = 0.0;
                    for (k = 0; k < 3; k++) {
                        dR[k] = periodicity(posc2[j][k] - posc1[i][k], Lbox[k]);
                        dRc += dR[k] * dR[k];
                    }
                    R = sqrt(dRc);
                    if (Cr < R && R < CR) {
                        for (m = 0; m < atoms1; m++) {
                            for (n = 0; n < atoms2; n++) {
                                drc = 0.0; r = 0.0;
                                for (g = 0; g < 3; g++) {
                                    dr[g] = periodicity(pos2[j][n][g] - pos1[i][m][g], Lbox[g]);
                                    drc += dr[g] * dr[g];
                                }
                                r = sqrt(drc);
                                VDW[m][n] = 4 * Etable[m][n] * (pow((Atable[m][n] / r), 12) - pow((Atable[m][n] / r), 6));
                                COU[m][n] = 138.9355 / D * atNCM[1][tN1][m] * atNCM[1][tN2][n] / r;
                                Vv0[i][j] += VDW[m][n];
                                Cc0[i][j] += COU[m][n];
                            }
                        }
                        Number[i] ++;
                        //if(i==1){printf("index=%5d Vvo= %8.4f Cc0= %8.4f Dis= %8.4f\n",index, Vv0[i][j],Cc0[i][j],R);}
                    }
                    else {
                        Vv0[i][j] = 0.0; Cc0[i][j] = 0.0;
                    }
                    Vv1[i] += Vv0[i][j];
                    Cc1[i] += Cc0[i][j];
                }

                NUM[meZ] += Number[i];

                if (Number[i] != 0) {
                    Vv2[meZ] += Vv1[i] / Number[i]; 
                    Cc2[meZ] += Cc1[i] / Number[i]; 
                }
                else {
                    Vv2[meZ] += Vv1[i];
                    Cc2[meZ] += Cc1[i];
                }
                //printf("\nTotal frame: %8.4e %8.4e\n",Vv2[meZ],Cc2[meZ]);
            }
        }
        for (i = 0; i < nbin; i++) {
            if (numbin[i] != 0) {
                V[index][i] += Vv2[i] / numbin[i];
                N[index][i] += double(NUM[i]) / numbin[i];
                C[index][i] += Cc2[i] / numbin[i];
                //printf("\nTotal frame: %8.4e %8.4e %8.4e\n",V[index][i],C[index][i],N[index][i]);		
            }
        }
        dt = tmpt;
        tmpt = fr.time;
        dt = tmpt - dt;
        index++;
        step++;
        if (step == endframe)
            break;
    } while (read_next_frame(oenv, status, &fr));

    sprintf(name0, "Energy-%d-%d-%.2f-%.2f-%.2f-%.2f-%d.dat", tN1 + 1, tN2 + 1, Cr, CR, lowPos, upPos, nbin);

    outputfileR = fopen(name0, "w");

    printf("Analysis done,4\n");

    for (i = 0; i < index; i++) { totV[i] = 0; totC[i] = 0; totN[i] = 0; }
    tV = 0; tN = 0; tC = 0;
    for (i = 0; i < nbin; i++) {
        for (j = 0; j < index; j++) {
            if (V[j][i] != 0) {
                totV[i] += V[j][i];
                totC[i] += C[j][i];
                totN[i] += N[j][i];
                T[i] ++;
            }
        }
    }
    printf("Analysis done,5\n");
    for (i = 0; i < nbin; i++) {
        tV += totV[i];
        tC += totC[i];
        tN += totN[i];
        fprintf(outputfileR, "%8.4e ", (i)*dbin);
        fprintf(outputfileR, "%8.4e %8.4e %8.4e %8.4e", totV[i] / T[i], totC[i] / T[i], totV[i] / T[i] + totC[i] / T[i], totN[i] / T[i]);
        fprintf(outputfileR, "\n");
    }
    printf("\nTotal frame:    %10d\n", index);
    printf("Total time(ps): %10.3f\n", fr.time - tmpt);
    printf("Analysis done\n");
    return 0;
}


//reference function
int load_alpha_epsilon(int atom1, int atom2, double **A, double **E, double *a1, double *a2, double *e1, double *e2, int LJtype)
{
    int m, n;
    if (LJtype == 1)
        for (m = 0; m < atom1; m++) {
            for (n = 0; n < atom2; n++) {
                A[m][n] = 0.5*(a1[m] + a2[n]);
                E[m][n] = sqrt(e1[m] * e2[n]);
            }
        }
    if (LJtype == 2)
        for (m = 0; m < atom1; m++) {
            for (n = 0; n < atom2; n++) {
                A[m][n] = sqrt((a1[m] * a2[n]));
                E[m][n] = sqrt(e1[m] * e2[n]);
            }
        }
    return 0;
}

int load_positionCenter(t_trxframe fr, int nm, int atoms, int *index, double **pos, double Lbox[3], double *deloc, double totNCM)
{
    int i, j, k, molN, m;
    double center[3], tmpPos, tmpCenter, halfLbox[3];
    for (k = 0; k < 3; k++)
        halfLbox[k] = Lbox[k] / 2;

    for (i = 0; i < nm; i++) {
        molN = index[i * atoms];

        for (k = 0; k < 3; k++) {
            tmpCenter = 0.0;  // initiate for each molecule
            for (j = 0; j < atoms; j++) {
                tmpPos = fr.x[molN + j][k];
                // first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
                while (tmpPos - fr.x[molN][k] > halfLbox[k])
                    tmpPos -= Lbox[k];
                while (tmpPos - fr.x[molN][k] < -halfLbox[k])
                    tmpPos += Lbox[k];

                tmpCenter += tmpPos * deloc[j];
            }
            pos[i][k] = tmpCenter / totNCM;
        }
    }

    return 0;
}

int load_position(t_trxframe fr, int nm, int atoms, int *index, double ***posmolecule, double Lbox[3])
{
    int i, j, k, molN;
    double center[3], tmpPos[3], halfLbox[3];
    for (k = 0; k < 3; k++)
        halfLbox[k] = Lbox[k] / 2;

    for (i = 0; i < nm; i++) {
        molN = index[i * atoms];
        for (j = 0; j < atoms; j++)
            for (k = 0; k < 3; k++) {
                tmpPos[k] = fr.x[molN + j][k];
                // first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
                while (tmpPos[k] - fr.x[molN][k] > halfLbox[k])
                    tmpPos[k] -= Lbox[k];
                while (tmpPos[k] - fr.x[molN][k] < -halfLbox[k])
                    tmpPos[k] += Lbox[k];

                posmolecule[i][j][k] = tmpPos[k];
            }
    }

    return 0;
}

double periodicity(double dx, double box)
{
    while (dx > box / 2.0)
        dx -= box;
    while (dx < -box / 2.0)
        dx += box;

    return dx;
}

double * CreateVector(int rows)
{
    return new double[rows]();
}


double ** CreateMatrix(int rows, int cols)
{
    int  i;
    double **m = new double*[rows];
    for (i = 0; i < rows; i++) {
        m[i] = new double[cols]();
    }
    return m;
}

int main(int argc, char *argv[])
{
    return gmx_run_cmain(argc, argv, &my_main);
}
