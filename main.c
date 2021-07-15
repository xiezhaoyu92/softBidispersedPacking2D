//
//  main.c
//  softBidispersity2D
//
//  Created by Zhaoyu Xie on 12/23/19.
//  Copyright © 2019 Zhaoyu Xie. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <stdlib.h>
#include "particles.h"
#include "operation.h"
#include "mt19937ar.h"

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON = 1e-15
#endif

int main(int argc, char * argv[]) {
    double L0 = 1;
    int np = 10;
    double rad = 0.1;
    double radRatio = 1; //larger to small
    double mixRatio = 0.5;
    
    double diffusionTimescale = 0.25;
    
    int animationQ = FALSE;
    int animationRate = 100;
    int readUnfinishedQ = FALSE;
    int energyRecordQ = TRUE;
    int diffusionQ = FALSE;
    
    static struct option longOptions[]={
        {"radius",              required_argument,  NULL, 'r'},
        {"radRatio",            required_argument,  NULL, 'R'},
        {"particleNumber",      required_argument,  NULL, 'n'},
        {"mixRatio",            required_argument,  NULL, 'i'},
        {"initialSize",         required_argument,  NULL, 'l'},
        {"animate",             no_argument,        NULL, 'm'},
        {"dtDiffusionScale",    required_argument,  NULL, 's'},
        {"readUnfinished",      no_argument,        NULL, 'u'},
        {"diffusion",           no_argument,        NULL, 'd'},
        {0,                     0,                  0,     0 }
    };
    int optIndex = 0;
    int opt;
    while ((opt = getopt_long(argc, argv,"r:R:n:i:l:ms:ud",
                              longOptions, &optIndex )) != -1) {
        switch (opt) {
            case 0:
                break;
            case 'r' :
                rad = atof(optarg);
                break;
            case 'R' :
                radRatio = atof(optarg);
                break;
            case 'n' :
                np = atoi(optarg);
                break;
            case 'i' :
                mixRatio = atof(optarg);
                break;
            case 'l' :
                L0 = atof(optarg);
                break;
            case 'm' :
                animationQ = TRUE;
                break;
            case 's':
                diffusionTimescale = atof(optarg);
                break;
            case 'u' :
                readUnfinishedQ = TRUE;
                break;
            case 'd' :
                diffusionQ = TRUE;
                break;
            default:
                exit(1);
        }
    }
    
    double L;
    double LOld;
    
    double simTimeStart = 0;
    double nextRelaxationTime = 0;
    int relaxationStep = 0;
    double simTime;
    int simStep = 0;
    
    if(readUnfinishedQ){
        FILE *radFile=NULL;
        radFile = fopen("rad.dat","r");
        if (radFile) {
            fscanf(radFile, "%lf", &rad);
            fclose(radFile);
        } else {
            printf("radFile pointer is null\n");
            exit(1);
        }
        FILE *npFile=NULL;
        npFile = fopen("npts.dat","r");
        if (npFile) {
            fscanf(npFile, "%i", &np);
            fclose(npFile);
        } else {
            printf("npFile pointer is null\n");
            exit(1);
        }
        FILE *rFileTemp=NULL;
        rFileTemp = fopen("sizeTemp.dat","r");
        if (rFileTemp) {
            fscanf(rFileTemp, "%lf", &L);
            fclose(rFileTemp);
        } else {
            printf("rFile pointer is null\n");
            exit(1);
        }
        FILE *rInFile=NULL;
        rInFile = fopen("sizeInit.dat", "r");
        if (rInFile) {
            fscanf(rInFile, "%lf", &L0);
            fclose(rInFile);
        } else {
            printf("rInFile pointer is null\n");
            exit(1);
        }
        FILE *simTimeTempFile=NULL;
        simTimeTempFile = fopen("simTimeTemp.dat","r");
        if (simTimeTempFile) {
            fscanf(simTimeTempFile, "%lf", &simTimeStart);
            fclose(simTimeTempFile);
        } else {
            printf("simTimeTempFile pointer is null\n");
            exit(1);
        }
        FILE *nextRelaxationFile=NULL;
        nextRelaxationFile = fopen("nextRelaxationTime.dat","r");
        if (nextRelaxationFile) {
            fscanf(nextRelaxationFile, "%lf", &nextRelaxationTime);
            fclose(nextRelaxationFile);
        } else {
            printf("nextRelaxationFile pointer is null\n");
            exit(1);
        }
        FILE *steps = NULL;
        steps = fopen("stepsTemp.dat","r");
        if(steps) {
            fscanf(steps, "%d", &relaxationStep);
            fclose(steps);
        }
        else {
            printf("steps pointer is null\n");
            exit(1);
        }
    }
    else {
        FILE *radFile=NULL;
        radFile = fopen("rad.dat","w");
        if (radFile) {
            fprintf(radFile, "%.15lf\n", rad);
            fclose(radFile);
        } else {
            printf("abFile pointer is null\n");
            exit(1);
        }
        FILE *npFile=NULL;
        npFile = fopen("npts.dat","w");
        if (npFile) {
            fprintf(npFile, "%i\n", np);
            fclose(npFile);
        } else {
            printf("npFile pointer is null\n");
            exit(1);
        }
        FILE *rInFile=NULL;
        rInFile = fopen("sizeInit.dat", "w");
        if (rInFile) {
            fprintf(rInFile, "%.15lf\n", L0);
            fclose(rInFile);
        } else {
            printf("rInFile pointer is null\n");
            exit(1);
        }
        
        L = L0;
        
    }
    
    double simTimeLastRelax = simTimeStart;
    double nextRelaxationTimeLastRelax = nextRelaxationTime;
    
    particle p[np];
    
    double dtDiffusion = 1;
    double D = pow(2*rad/1000, 2)/2; //diffusion coefficient along long axis
    double delta = sqrt(2*D);
    double simTimeFinal = 2*rad*rad/(D*diffusionTimescale);
    double relaxationSteps = 1e5;
    double dtRelaxationMax = simTimeFinal/relaxationSteps;
    double dtRelaxation = dtRelaxationMax;
    double dtTol = 1e-5*dtRelaxation;
    int outputSteps = 1;
/*    if (dtDiffusion>dtRelaxation) {
        printf("Fast relaxation rate: using smaller timestep.\n");
        dtDiffusion = dtRelaxation;
        //printf("%e\n",dtDiffusion0);
    }
*/
    if(readUnfinishedQ) {
        FILE *dtTempFile=NULL;
        dtTempFile = fopen("dtTemp.dat","r");
        if (dtTempFile) {
            fscanf(dtTempFile, "%lf", &dtRelaxation);
            fclose(dtTempFile);
        } else {
            printf("dtTempFile pointer is null\n");
            exit(1);
        }
    }
    nextRelaxationTime += dtRelaxation;
    
    double configEnergy = 0;
    
    int nPart[2];
    nPart[0] = (int)(L/(2*rad));
    nPart[1] = (int)(L/(2*rad));
    double partBoundX[nPart[0]];
    double partBoundY[nPart[1]];
    for(int i=0; i<nPart[0]; i++)
        partBoundX[i] = i*2*rad;
    for(int i=0; i<nPart[1]; i++)
        partBoundY[i] = i*2*rad;

    clock_t begin, end;
    double timeSpent;
    randinitialize();
    begin = clock();
    
    if(readUnfinishedQ){
        FILE *configurationTemp=NULL;
        FILE *radiusFile=NULL;
        configurationTemp = fopen("configurationTemp.asc","r");
        radiusFile = fopen("radii.dat", "r");
        if (configurationTemp&&radiusFile) {
            for (int j = 0; j < np; j++) {
                fscanf(configurationTemp, "%lf %lf", &(p[j].position[0]), &(p[j].position[1]));
                fscanf(radiusFile, "%lf", &(p[j].rad));
                partitionOneParticle(&(p[j]), partBoundX, partBoundY, nPart);
                for(int k=0; k<2; k++) {
                    p[j].postPreRelaxPosition[k] = p[j].position[k];
                }
                for(int k=0; k<2; k++) {
                    p[j].postPreRelaxCoord[k] = p[j].coord[k];
                }
            }
            fclose(configurationTemp);
        }
        else {
            printf("configuration or radii pointer is null\n");
            exit(1);
        }
    }
    else {
        for (int i=0; i< (int)(np*mixRatio); i++) {
            p[i].rad=rad;
        }
        for (int i=(int)(np*mixRatio); i<np; i++) {
            p[i].rad=rad/radRatio;
        }
        FILE *radiusFile=NULL;
        radiusFile = fopen("radii.dat","w");
        if (radiusFile) {
            for(int i=0; i<np; i++)
                fprintf(radiusFile, "%lf\n", p[i].rad);
            fclose(radiusFile);
        } else {
            printf("radiusFile pointer is null\n");
            exit(1);
        }
        FILE *init=NULL; // after projection
        init = fopen("init.asc","w");
        if (!init) {
            printf("file init.asc failed to open");
            exit(1);
        }
        for(int i=0; i<np; i++) {
            int overlapCheck;
            do {
                overlapCheck = FALSE;
                double x = L*genrand_real2();
                double y = L*genrand_real2();
                p[i].position[0] = x;
                p[i].position[1] = y;
                partitionOneParticle(&(p[i]), partBoundX, partBoundY, nPart);
                for(int j=0; j<i; j++) {
                    overlapCheck = overlapQ(&p[i], &p[j], nPart, L);
                    if(overlapCheck)
                        break;
                }
            }while(overlapCheck);
            fprintf(init, "%.15lf %.15lf\n", p[i].position[0], p[i].position[1]);
        }
        if(init) fclose(init);
    }
    
    FILE *animationFile=NULL;
    if (animationQ) {
        animationFile = fopen("animation.dat", "a");
        if (animationFile) {
            fprintf(animationFile, "%.15lf\n", L);
            for (int j=0; j<np; j++)
                fprintf(animationFile, "%.15lf %.15lf\n", p[j].position[0], p[j].position[1]);
            fclose(animationFile);
        }
        else {
            printf("animationFile pointer is null\n");
            exit(1);
        }
    }
    
    double Vmin = 1e-16;
    
    FILE *energyFile = NULL;
    FILE *configFile = NULL;
    
    configEnergy = totalRepulsiveSpringEnergy(p, np, L, nPart);
    
    int rollbackQ = FALSE;
    double dtRelaxationOld = 0;
    for (simTime = simTimeStart; simTime <= simTimeFinal; simTime += dtRelaxation) {
        if (simTime >= nextRelaxationTime) {
            LOld = L;
            L = L0*(1-simTime/simTimeFinal);
            updatePartitions(nPart, partBoundX, partBoundY, L, rad);
            for(int j=0; j<np; j++) {
                projectIntoNewArea(&(p[j]), L, LOld);
                partitionOneParticle(&(p[j]), partBoundX, partBoundY, nPart);
            }
            
            for(int j=0; j<np; j++)
                for(int k=0; k<2; k++)
                    p[j].force[k] = 0;
            //conjugateGradientDescentGoldenSearch(p, np, rad, L, partBoundX, partBoundY, nPart);
            conjugateGradientDescentGoldenSearchSteps(p, np, rad, L, partBoundX, partBoundY, nPart, 1e6);
            
            configEnergy = totalRepulsiveSpringEnergy(p, np, L, nPart);
/*
            if(configEnergy/np>=Vmin) {
                double LL = L*0.99;
                updatePartitions(nPart, partBoundX, partBoundY, LL, rad);
                for(int i=0; i<np; i++) {
                    projectIntoNewArea(&(p[i]), LL, L);
                    partitionOneParticle(&(p[i]), partBoundX, partBoundY, nPart);
                }
                for(int j=0; j<np; j++)
                    for(int k=0; k<2; k++)
                        p[j].force[k] = 0;
                conjugateGradientDescentGoldenSearch(p, np, rad, LL, partBoundX, partBoundY, nPart);
                updatePartitions(nPart, partBoundX, partBoundY, L, rad);
                for(int i=0; i<np; i++) {
                    projectIntoNewArea(&(p[i]), L, LL);
                    partitionOneParticle(&(p[i]), partBoundX, partBoundY, nPart);
                }
                for(int j=0; j<np; j++)
                    for(int k=0; k<2; k++)
                        p[j].force[k] = 0;
                conjugateGradientDescentGoldenSearch(p, np, rad, L, partBoundX, partBoundY, nPart);
                configEnergy = totalRepulsiveSpringEnergy(p, np, L, nPart);
            }
*/            
            if(energyRecordQ){
                energyFile = fopen("energy.dat","a");
                if(energyFile) {
                    fprintf(energyFile, "%.15lf %.20lf\n", simTime, configEnergy);
                    fclose(energyFile);
                }
                else {
                    printf("energyFile pointer is null\n");
                    exit(1);
                }
            }
            
            if(configEnergy/np<Vmin&&LOld!=L) {
                for (int j=0; j<np; j++) {
                    for(int k=0; k<2; k++) {
                        p[j].postPreRelaxPosition[k] = p[j].position[k];
                    }
                    for(int k=0; k<2; k++) {
                        p[j].postPreRelaxCoord[k] = p[j].coord[k];
                    }
                }
                simTimeLastRelax = simTime;
                nextRelaxationTimeLastRelax = nextRelaxationTime;
                relaxationStep++;
                
                if(rollbackQ) {
                    dtRelaxation = dtRelaxation/2;
                    rollbackQ = FALSE;
                    dtRelaxationOld = 0;
                }
                else {
                    if(dtRelaxation == dtRelaxationOld)
                        dtRelaxation = dtRelaxation*2;
                    else
                        dtRelaxationOld = dtRelaxation;
                }
                if(dtRelaxation>dtRelaxationMax) {
                    dtRelaxation = dtRelaxationMax;
                    dtRelaxationOld = dtRelaxation;
                }
                
                if (animationQ && relaxationStep % animationRate == 0) {
                    animationFile = fopen("animation.dat", "a");
                    if (animationFile) {
                        fprintf(animationFile, "%.15lf\n", L);
                        for (int j=0; j<np; j++)
                            fprintf(animationFile, "%.15lf %.15lf\n", p[j].position[0], p[j].position[1]);
                        fclose(animationFile);
                    } else
                        printf("animationFile pointer is null\n");
                }
                
                if (relaxationStep%outputSteps==0) {
                    FILE *configurationTemp=NULL;
                    configurationTemp = fopen("configurationTemp.asc.tmp","w");
                    if (configurationTemp) {
                        for (int j = 0; j < np; j++)
                            fprintf(configurationTemp, "%.15lf %.15lf\n", p[j].position[0], p[j].position[1]);
                        fclose(configurationTemp);
                    }
                    else {
                        printf("configuration pointer is null\n");
                        exit(1);
                    }
                    FILE *rFileTemp=NULL;
                    rFileTemp = fopen("sizeTemp.dat.tmp","w");
                    if (rFileTemp) {
                        fprintf(rFileTemp, "%.15lf\n", L);
                        fclose(rFileTemp);
                    } else {
                        printf("rFile pointer is null\n");
                        exit(1);
                    }
                    FILE *simTimeTempFile=NULL;
                    simTimeTempFile = fopen("simTimeTemp.dat.tmp","w");
                    if (simTimeTempFile) {
                        fprintf(simTimeTempFile, "%.15lf\n", simTime);
                        fclose(simTimeTempFile);
                    } else {
                        printf("simTimeTempFile pointer is null\n");
                        exit(1);
                    }
                    FILE *nextRelaxationFile=NULL;
                    nextRelaxationFile = fopen("nextRelaxationTime.dat.tmp","w");
                    if (nextRelaxationFile) {
                        fprintf(nextRelaxationFile, "%.15lf\n", nextRelaxationTime);
                        fclose(nextRelaxationFile);
                    } else {
                        printf("nextRelaxationFile pointer is null\n");
                        exit(1);
                    }
                    FILE *stepsTemp = NULL;
                    stepsTemp = fopen("steps.dat.tmp","w");
                    if(stepsTemp) {
                        fprintf(stepsTemp, "%d\n", relaxationStep);
                        fclose(stepsTemp);
                    }
                    else {
                        printf("stepsTemp pointer is null\n");
                        exit(1);
                    }
                    
                    system("mv -f configurationTemp.asc.tmp configurationTemp.asc");
                    system("mv -f sizeTemp.dat.tmp sizeTemp.dat");
                    system("mv -f simTimeTemp.dat.tmp simTimeTemp.dat");
                    system("mv -f nextRelaxationTime.dat.tmp nextRelaxationTime.dat");
                    system("mv -f steps.dat.tmp stepsTemp.dat");
                    //system("mv -f dtTemp.dat.tmp dtTemp.dat");
                }
                
            }
            else {
                if(configEnergy/np<2*Vmin)
                    break;
                rollbackQ = TRUE;
                dtRelaxation = dtRelaxation/2;
                //if (dtDiffusion>dtRelaxation)
                //    dtDiffusion = dtRelaxation;
                L = LOld;
                for (int j=0; j<np; j++) {
                    for(int k=0; k<2; k++) {
                        p[j].position[k] = p[j].postPreRelaxPosition[k];
                    }
                    for(int k=0; k<2; k++) {
                        p[j].coord[k] = p[j].postPreRelaxCoord[k];
                    }
                }
                updatePartitions(nPart, partBoundX, partBoundY, L, rad);
                simTime = simTimeLastRelax;
                nextRelaxationTime = nextRelaxationTimeLastRelax;
            }
            if (relaxationStep%outputSteps==0) {
                FILE *dtTempFile=NULL;
                dtTempFile = fopen("dtTemp.dat","w");
                if (dtTempFile) {
                    fprintf(dtTempFile, "%.15lf\n", dtRelaxation);
                    fclose(dtTempFile);
                } else {
                    printf("dtTempFile pointer is null\n");
                    exit(1);
                }
            }
            nextRelaxationTime += dtRelaxation;
        }
        
        simStep++;
        printf("simTime/simTimeFinal: %.16lf\n", simTime/simTimeFinal);
    }
    
    if (animationQ) {
        animationFile = fopen("animation.dat", "a");
        if (animationFile) {
            fprintf(animationFile, "%.15lf\n", L);
            for (int j=0; j<np; j++)
                fprintf(animationFile, "%.15lf %.15lf\n", p[j].position[0], p[j].position[1]);
            fclose(animationFile);
        } else
            printf("animationFile pointer is null\n");
    }
    
    end = clock();
    timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;

    FILE *timeFile=NULL;
    timeFile = fopen("timeS.dat","w");
    if (timeFile) {
        fprintf(timeFile, "%lf\n", timeSpent);
        fclose(timeFile);
    } else {
        printf("timeFile pointer is null\n");
        exit(1);
    }
    
    FILE *rFile=NULL;
    rFile = fopen("sizeS.dat","w");
    if (rFile) {
        fprintf(rFile, "%.15lf\n", L);
        fclose(rFile);
    } else {
        printf("rFile pointer is null\n");
        exit(1);
    }
    
    FILE *configuration=NULL;
    configuration = fopen("configurationS.asc","w");
    if (configuration) {
        for (int j = 0; j < np; j++)
            fprintf(configuration, "%.15lf %.15lf\n", p[j].position[0], p[j].position[1]);
        fclose(configuration);
    }
    else {
        printf("configuration pointer is null\n");
        exit(1);
    }
    
    FILE *simTimeArrestedFile=NULL;
    simTimeArrestedFile = fopen("simTimeArrestedS.dat","w");
    if (simTimeArrestedFile) {
        fprintf(simTimeArrestedFile, "%.15lf\n", simTime);
        fclose(simTimeArrestedFile);
    }
    else {
        printf("simTimeArrestedFile pointer is null\n");
        exit(1);
    }
    
    FILE *nextRelaxationArrestedFile=NULL;
    nextRelaxationArrestedFile = fopen("nextRelaxationTimeArrestedS.dat","w");
    if (nextRelaxationArrestedFile) {
        fprintf(nextRelaxationArrestedFile, "%.15lf\n", nextRelaxationTime);
        fclose(nextRelaxationArrestedFile);
    }
    else {
        printf("nextRelaxationArrestedFile pointer is null\n");
        exit(1);
    }
    
    FILE *steps = NULL;
    steps = fopen("stepsS.dat","w");
    if(steps) {
        fprintf(steps, "%d\n", relaxationStep);
        fclose(steps);
    }
    else {
        printf("steps pointer is null\n");
        exit(1);
    }
    
    printf("Hello, World!\n");
    return 0;
}
