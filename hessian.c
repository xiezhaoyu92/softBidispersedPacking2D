//
//  hessian.c
//  softEllipsoid2D
//
//  Created by Zhaoyu Xie on 11/19/19.
//  Copyright © 2019 Zhaoyu Xie. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hessian.h"
#include "operation.h"
#include "particles.h"
#include "solveHessian.h"

double Vii(int i, particle p1, particle p2) {
    double epsilon = 1;
    double r[2];
    vectorSubtract(p1.position, p2.position, r);
    double R = vectorNorm(r);
    double sigma = p1.rad+p2.rad;
    return epsilon/(sigma*sigma)*pow(p1.position[i]-p2.position[i],2)/(R*R)-epsilon/sigma*(1-R/sigma)*(1/R-pow(p1.position[i]-p2.position[i],2)/(R*R*R));
}
double Vij(int i, int j, particle p1, particle p2) {
    double epsilon = 1;
    double r[2];
    vectorSubtract(p1.position, p2.position, r);
    double R = vectorNorm(r);
    double sigma = p1.rad+p2.rad;
    return epsilon/(sigma*sigma)*(p1.position[i]-p2.position[i])*(p1.position[j]-p2.position[j])/(R*R)+epsilon/sigma*(1-R/sigma)*(p1.position[i]-p2.position[i])*(p1.position[j]-p2.position[j])/(R*R*R);
}
double ViVi(int i, particle p1, particle p2) {
    double epsilon = 1;
    double r[2];
    vectorSubtract(p1.position, p2.position, r);
    double R = vectorNorm(r);
    double sigma = p1.rad+p2.rad;
    return -epsilon/(sigma*sigma)*pow(p1.position[i]-p2.position[i],2)/(R*R)-epsilon/sigma*(1-R/sigma)*(-1/R+pow(p1.position[i]-p2.position[i],2)/(R*R*R));
}
double ViVj(int i, int j, particle p1, particle p2) {
    double epsilon = 1;
    double r[2];
    vectorSubtract(p1.position, p2.position, r);
    double R = vectorNorm(r);
    double sigma = p1.rad+p2.rad;
    return -epsilon/(sigma*sigma)*(p1.position[i]-p2.position[i])*(p1.position[j]-p2.position[j])/(R*R)-epsilon/sigma*(1-R/sigma)*(p1.position[i]-p2.position[i])*(p1.position[j]-p2.position[j])/(R*R*R);
}

double secondTotalD(int m, int i, particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[]) {
    double result = 0;
    particle p1 = p[m];
    for(int n=0; n<np; n++) {
        if(n==m||!overlapQ(&(p[m]), &(p[n]), nPart, L))
            continue;
        
        particle p2 = p[n];
        if(p2.coord[0]-p1.coord[0]==nPart[0]-1)
            p2.position[0] = p2.position[0]-L;
        else if(p2.coord[0]-p1.coord[0]==-(nPart[0]-1))
            p2.position[0] = p2.position[0]+L;
        
        if(p2.coord[1]-p1.coord[1]==nPart[1]-1)
            p2.position[1] = p2.position[1]-L;
        else if(p2.coord[1]-p1.coord[1]==-(nPart[1]-1))
            p2.position[1] = p2.position[1]+L;
        
        result += Vii(i, p1, p2);
    }
    return result;
}

double secondPartialDSingle(int m, int i, int j, particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[]){
    double result = 0;
    particle p1 = p[m];
    for(int n=0; n<np; n++) {
        if(n==m||!overlapQ(&(p[m]), &(p[n]), nPart, L))
            continue;
        
        particle p2 = p[n];
        if(p2.coord[0]-p1.coord[0]==nPart[0]-1)
            p2.position[0] = p2.position[0]-L;
        else if(p2.coord[0]-p1.coord[0]==-(nPart[0]-1))
            p2.position[0] = p2.position[0]+L;
        
        if(p2.coord[1]-p1.coord[1]==nPart[1]-1)
            p2.position[1] = p2.position[1]-L;
        else if(p2.coord[1]-p1.coord[1]==-(nPart[1]-1))
            p2.position[1] = p2.position[1]+L;
        
        result += Vij(i, j, p1, p2);
    }
    return result;
}

double secondPartialDPair(int m, int n, int i, int j, particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[]) {
    particle p1 = p[m];
    particle p2 = p[n];
    if(p2.coord[0]-p1.coord[0]==nPart[0]-1)
        p2.position[0] = p2.position[0]-L;
    else if(p2.coord[0]-p1.coord[0]==-(nPart[0]-1))
        p2.position[0] = p2.position[0]+L;
    
    if(p2.coord[1]-p1.coord[1]==nPart[1]-1)
        p2.position[1] = p2.position[1]-L;
    else if(p2.coord[1]-p1.coord[1]==-(nPart[1]-1))
        p2.position[1] = p2.position[1]+L;
    
    if(i==j)
        return ViVi(i, p1, p2);
    else
        return ViVj(i, j, p1, p2);
}

double repulsiveSpringEnergyOneParticle(int i, particle pp, particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[]) {
    
    double energy = 0;
    for(int j=0; j<np; j++){
        if(j==i||!overlapQ(&(p[i]), &(p[j]), nPart, L))
            continue;
        double x2;
        double y2;
        
        if(p[j].coord[0]==pp.coord[0]||abs(p[j].coord[0]-pp.coord[0])==1)
            x2 = p[j].position[0];
        else if(p[j].coord[0]-pp.coord[0]==nPart[0]-1)
            x2 = p[j].position[0]-L;
        else if(p[j].coord[0]-pp.coord[0]==-(nPart[0]-1))
            x2 = p[j].position[0]+L;
        else
            continue;
        
        if(p[j].coord[1]==pp.coord[1]||abs(p[j].coord[1]-pp.coord[1])==1)
            y2 = p[j].position[1];
        else if(p[j].coord[1]-pp.coord[1]==nPart[1]-1)
            y2 = p[j].position[1]-L;
        else if(p[j].coord[1]-pp.coord[1]==-(nPart[1]-1))
            y2 = p[j].position[1]+L;
        else
            continue;
        
        energy = energy+repulsiveSpringPotential(pp.position[0], pp.position[1], x2, y2, pp.rad, p[j].rad);
    }
    return energy;
}

double repulsiveSpringEnergyPairs(particle p1, particle p2, int np, double L, int nPart[], double partBoundX[], double partBoundY[]) {
    double x1 = p1.position[0];
    double y1 = p1.position[1];
    double x2;
    double y2;
    
    if(p2.coord[0]==p1.coord[0]||abs(p2.coord[0]-p1.coord[0])==1)
        x2 = p2.position[0];
    else if(p2.coord[0]-p1.coord[0]==nPart[0]-1)
        x2 = p2.position[0]-L;
    else if(p2.coord[0]-p1.coord[0]==-(nPart[0]-1))
        x2 = p2.position[0]+L;
    else
        return 0;
    
    if(p2.coord[1]==p1.coord[1]||abs(p2.coord[1]-p1.coord[1])==1)
        y2 = p2.position[1];
    else if(p2.coord[1]-p1.coord[1]==nPart[1]-1)
        y2 = p2.position[1]-L;
    else if(p2.coord[1]-p1.coord[1]==-(nPart[1]-1))
        y2 = p2.position[1]+L;
    else
        return 0;
    
    return repulsiveSpringPotential(x1, y1, x2, y2, p1.rad, p2.rad);
}
/*
double secondTotalD(int m, int i, particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[]) {
    double delta = 0.000001;
    particle p0=p[m];
    particle p1=p[m];
    particle p2=p[m];
    
    p1.position[i] += delta;
    p2.position[i] -= delta;
    while(p1.position[i]>L)
        p1.position[i] = p1.position[i]-L;
    while(p2.position[i]<0)
        p2.position[i] = p2.position[i]+L;
    partitionOneParticle(&p1, partBoundX, partBoundY, nPart);
    partitionOneParticle(&p2, partBoundX, partBoundY, nPart);
    
    double E0 = repulsiveSpringEnergyOneParticle(m, p0, p, np, L, nPart, partBoundX, partBoundY);
    double E1 = repulsiveSpringEnergyOneParticle(m, p1, p, np, L, nPart, partBoundX, partBoundY);
    double E2 = repulsiveSpringEnergyOneParticle(m, p2, p, np, L, nPart, partBoundX, partBoundY);
    return (E1+E2-2*E0)/(delta*delta);
}

//j is larger than i, since only j may refer to angle
double secondPartialDSingle(int m, int i, int j, particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[]) {
    double delta = 0.000001;
    particle p1=p[m];
    particle p2=p[m];
    particle p3=p[m];
    particle p4=p[m];
    
    p1.position[i] += delta;
    while(p1.position[i]>L)
        p1.position[i] = p1.position[i]-L;
    p1.position[j] += delta;
    while(p1.position[j]>L)
        p1.position[j] = p1.position[j]-L;
    partitionOneParticle(&p1, partBoundX, partBoundY, nPart);
        
    p2.position[i] -= delta;
    while(p2.position[i]<0)
        p2.position[i] = p1.position[i]+L;
    p2.position[j] -= delta;
    while(p2.position[j]<0)
        p2.position[j] = p2.position[j]+L;
    partitionOneParticle(&p2, partBoundX, partBoundY, nPart);
        
    p3.position[i] += delta;
    while(p3.position[i]>L)
        p3.position[i] = p3.position[i]-L;
    p3.position[j] -= delta;
    while(p3.position[j]<0)
        p3.position[j] = p3.position[j]+L;
    partitionOneParticle(&p3, partBoundX, partBoundY, nPart);
        
    p4.position[i] -= delta;
    while(p4.position[i]<0)
        p4.position[i] = p4.position[i]+L;
    p4.position[j] += delta;
    while(p4.position[j]>L)
        p4.position[j] = p4.position[j]-L;
    partitionOneParticle(&p4, partBoundX, partBoundY, nPart);
    
    double E1 = repulsiveSpringEnergyOneParticle(m, p1, p, np, L, nPart, partBoundX, partBoundY);
    double E2 = repulsiveSpringEnergyOneParticle(m, p2, p, np, L, nPart, partBoundX, partBoundY);
    double E3 = repulsiveSpringEnergyOneParticle(m, p3, p, np, L, nPart, partBoundX, partBoundY);
    double E4 = repulsiveSpringEnergyOneParticle(m, p4, p, np, L, nPart, partBoundX, partBoundY);
    
    return (E1+E2-E3-E4)/(4*delta*delta);
}

double secondPartialDPair(int m, int n, int i, int j, particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[]) {
    double delta = 0.000001;
    particle pm1 = p[m];
    particle pm2 = p[m];
    particle pn1 = p[n];
    particle pn2 = p[n];
    
    pm1.position[i] += delta;
    pm2.position[i] -= delta;
    while(pm1.position[i]>L)
        pm1.position[i] = pm1.position[i]-L;
    while(pm2.position[i]<0)
        pm2.position[i] = pm2.position[i]+L;
    partitionOneParticle(&pm1, partBoundX, partBoundY, nPart);
    partitionOneParticle(&pm2, partBoundX, partBoundY, nPart);
    
    pn1.position[j] += delta;
    pn2.position[j] -= delta;
    while(pn1.position[j]>L)
        pn1.position[j] = pn1.position[j]-L;
    while(pn2.position[j]<0)
        pn2.position[j] = pn2.position[j]+L;
    partitionOneParticle(&pn1, partBoundX, partBoundY, nPart);
    partitionOneParticle(&pn2, partBoundX, partBoundY, nPart);
    
    double E1 = repulsiveSpringEnergyPairs(pm1, pn1, np, L, nPart, partBoundX, partBoundY);
    double E2 = repulsiveSpringEnergyPairs(pm2, pn2, np, L, nPart, partBoundX, partBoundY);
    double E3 = repulsiveSpringEnergyPairs(pm1, pn2, np, L, nPart, partBoundX, partBoundY);
    double E4 = repulsiveSpringEnergyPairs(pm2, pn1, np, L, nPart, partBoundX, partBoundY);
    
    return (E1+E2-E3-E4)/(4*delta*delta);
}
*/
int indexColumn(int np, int row, int column) {
    return column*2*np+row;
}
void hessianMatrixSphere(particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[], double *hessian) {
    for(int m=0; m<np; m++) {
        hessian[indexColumn(np,2*m,2*m)] = secondTotalD(m, 0, p, np, L, nPart, partBoundX, partBoundY);
        hessian[indexColumn(np,2*m+1,2*m+1)] = secondTotalD(m, 1, p, np, L, nPart, partBoundX, partBoundY);
        hessian[indexColumn(np,2*m,2*m+1)] = secondPartialDSingle(m, 0, 1, p, np, L, nPart, partBoundX, partBoundY);
        hessian[indexColumn(np,2*m+1,2*m)] = hessian[indexColumn(np,2*m,2*m+1)];
        for(int n=m+1; n<np; n++) {
            if(!overlapQ(&(p[m]), &(p[n]), nPart, L))
                continue;
            hessian[indexColumn(np,2*m,2*n)] = secondPartialDPair(m, n, 0, 0, p, np, L, nPart, partBoundX, partBoundY);
            hessian[indexColumn(np,2*m,2*n+1)] = secondPartialDPair(m, n, 0, 1, p, np, L, nPart, partBoundX, partBoundY);
            hessian[indexColumn(np,2*m+1,2*n)] = secondPartialDPair(m, n, 1, 0, p, np, L, nPart, partBoundX, partBoundY);
            hessian[indexColumn(np,2*m+1,2*n+1)] = secondPartialDPair(m, n, 1, 1, p, np, L, nPart, partBoundX, partBoundY);
            hessian[indexColumn(np,2*n,2*m)] = hessian[indexColumn(np,2*m,2*n)];
            hessian[indexColumn(np,2*n+1,2*m)] = hessian[indexColumn(np,2*m,2*n+1)];
            hessian[indexColumn(np,2*n,2*m+1)] = hessian[indexColumn(np,2*m+1,2*n)];
            hessian[indexColumn(np,2*n+1,2*m+1)] = hessian[indexColumn(np,2*m+1,2*n+1)];
        }
    }
}

void solveHessianSphere(particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[], int *info) {
    int nn=0;
    int pi[np];
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++)
            if(overlapQ(&(p[i]), &(p[j]), nPart, L)){
                if(nn==0){
                    pi[0]=i;
                    pi[1]=j;
                    nn=2;
                    continue;
                }
                int k;
                for(k=0; k<nn; k++)
                    if(pi[k]==i)
                        break;
                if(k==nn){
                    pi[k]=i;
                    nn++;
                }
                for(k=0; k<nn; k++)
                    if(pi[k]==j)
                        break;
                if(k==nn){
                    pi[k]=j;
                    nn++;
                }
            }
    particle pp[nn];
    for(int i=0; i<nn; i++)
        pp[i] = p[pi[i]];
    
    double *hessian = malloc(2*nn*2*nn*sizeof(double));
    for(int i=0; i<2*nn*2*nn; i++)
        hessian[i] = 0;
    hessianMatrixSphere(pp, nn, L, nPart, partBoundX, partBoundY, hessian);
    /*FILE *hes=fopen("hessian.txt","w");
    for(int i=0; i<2*nn*2*nn; i++)
        fprintf(hes,"%lf\n",hessian[i]);
    fclose(hes);*/
    
    for(int j=0; j<np; j++)
        for(int k=0; k<2; k++)
            p[j].force[k] = 0;
    addRepulsiveSpringForce(p, np, L, nPart);
    double forces[2*nn];
    for(int i=0; i<nn; i++) {
        forces[2*i] = p[pi[i]].force[0];
        forces[2*i+1] = p[pi[i]].force[1];
    }
    /*FILE *force=fopen("force.txt","w");
    for(int i=0; i<2*nn; i++)
        fprintf(force,"%.15lf\n",forces[i]);
    fclose(force);*/
    
    /*double oldForces[2*nn];
    for(int i=0; i<2*nn; i++)
        oldForces[i]=forces[i];
    double oldHessian[2*nn*2*nn];
    for(int i=0; i<2*nn*2*nn; i++)
        oldHessian[i]=hessian[i];*/
    
    int n = 2*nn;
    int nrhs = 1;
    int lda = n;
    int ldb = n;
    int ipiv[n];
    dgesv_(&n,&nrhs,hessian,&lda,ipiv,forces,&ldb,info);
    /*for(int i=0; i<2*nn; i++)
        printf("%d: %.15lf\n", i, forces[i]);*/

    /*for(int i=0; i<2*nn; i++){
        double sum=0;
        for(int j=0; j<2*nn; j++)
            sum += oldHessian[i*2*nn+j]*forces[j];
        printf("%.10lf, %.10lf, %.10lf\n",sum,oldForces[i],sum-oldForces[i]);
    }*/
    
    for(int j=0; j<np; j++)
        for(int k=0; k<4; k++)
            p[j].force[k] = 0;
    for(int i=0; i<nn; i++){
        p[pi[i]].force[0]=forces[2*i];
        p[pi[i]].force[1]=forces[2*i+1];
    }
    free(hessian);
}



/*
int  main() {
    
    int np;
    double L, rad;
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
    rFileTemp = fopen("sizeS.dat","r");
    if (rFileTemp) {
        fscanf(rFileTemp, "%lf", &L);
        fclose(rFileTemp);
    } else {
        printf("rFile pointer is null\n");
        exit(1);
    }
    FILE *radFile=NULL;
    radFile = fopen("rad.dat","r");
    if (radFile) {
        fscanf(radFile, "%lf", &rad);
        fclose(radFile);
    } else {
        printf("radFile pointer is null\n");
        exit(1);
    }
    
    particle p[np];
    int nPart[2];
    nPart[0] = (int)(L/(2*rad));
    nPart[1] = (int)(L/(2*rad));
    double partBoundX[nPart[0]];
    double partBoundY[nPart[1]];
    for(int i=0; i<nPart[0]; i++)
        partBoundX[i] = i*2*rad;
    for(int i=0; i<nPart[1]; i++)
        partBoundY[i] = i*2*rad;
    
    FILE *configurationTemp=NULL;
    FILE *radiusFile=NULL;
    configurationTemp = fopen("configurationS.asc","r");
    radiusFile = fopen("radii.dat", "r");
    if (configurationTemp&&radiusFile) {
        for (int j = 0; j < np; j++) {
            fscanf(configurationTemp, "%lf %lf", &(p[j].position[0]), &(p[j].position[1]));
            fscanf(radiusFile, "%lf", &(p[j].rad));
            partitionOneParticle(&(p[j]), partBoundX, partBoundY, nPart);
        }
        fclose(configurationTemp);
    }
    else {
        printf("configuration or radii pointer is null\n");
        exit(1);
    }
    
    double *hessian = malloc(2*np*2*np*sizeof(double));
    for(int i=0; i<2*np*2*np; i++)
        hessian[i] = 0;
    hessianMatrixSphere(p, np, L, nPart, partBoundX, partBoundY, hessian);
    for(int i=0;i<2*np;i++)
        printf("%d: %lf\n", i, hessian[i]);
    
    free(hessian);
    
    int info;
   
    solveHessianSphere(p, np, L, nPart, partBoundX, partBoundY, &info);
    printf("info: %d\n", info);
    FILE *forceFile = fopen("displacement.txt","w");
    for(int i=0;i<np;i++)
        fprintf(forceFile, "%.15lf %.15lf %.15lf\n", p[i].force[0],p[i].force[1]);
    fclose(forceFile);
}
*/
