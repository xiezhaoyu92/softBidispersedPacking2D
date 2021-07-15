//
//  particles.c
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#include "particles.h"
#include <math.h>
#include "operation.h"
#include <stdlib.h>
#include "hessian.h"

int overlapQ(particle *p1, particle *p2, int nPart[], double L) {
    particle temp;
    temp.rad = p2->rad;
    if(p2->coord[0]==p1->coord[0]||abs(p2->coord[0]-p1->coord[0])==1)
        temp.position[0]=p2->position[0];
    else if(p2->coord[0]-p1->coord[0]==nPart[0]-1)
        temp.position[0]=p2->position[0]-L;
    else if(p2->coord[0]-p1->coord[0]==-(nPart[0]-1))
        temp.position[0]=p2->position[0]+L;
    else
        return FALSE;
    if(p2->coord[1]==p1->coord[1]||abs(p2->coord[1]-p1->coord[1])==1)
        temp.position[1]=p2->position[1];
    else if(p2->coord[1]-p1->coord[1]==nPart[1]-1)
        temp.position[1]=p2->position[1]-L;
    else if(p2->coord[1]-p1->coord[1]==-(nPart[1]-1))
        temp.position[1]=p2->position[1]+L;
    else
        return FALSE;
    
    double r[2];
    vectorSubtract(p1->position, temp.position, r);
    if(vectorNorm(r)<(p1->rad+temp.rad))
        return TRUE;
    else
        return FALSE;
}

//move particle, and project onto the surface. Update the information of particle
void gradDescentConstrainedIntegrate(particle *p, double dt, double L) {
    //move particle
    for(int i=0; i<2; i++)
        p->position[i] = p->position[i] + (p->force[i])*dt + (p->forceStoc[i])*sqrt(dt);
    while(p->position[0]<0)
        p->position[0] = p->position[0]+L;
    while(p->position[0]>L)
        p->position[0] = p->position[0]-L;
    while(p->position[1]<0)
        p->position[1] = p->position[1]+L;
    while(p->position[1]>L)
        p->position[1] = p->position[1]-L;
}

//Add stochastic force
//rad is the radius of larger sphere
void addStochasticForce(double D, double rad, particle *p) {
    for(int i=0; i<2; i++)
        p->forceStoc[i] = sqrt(2*D)/rad*(p->rad)*randNormal();
}

void partitionOneParticle(particle *p, double partBoundX[], double partBoundY[], int nPart[]) {
    int i=0, j=0;
    while(i<nPart[0]&&p->position[0]>=partBoundX[i])
        i++;
    p->coord[0] = i-1;
    while(j<nPart[1]&&p->position[1]>=partBoundY[j])
        j++;
    p->coord[1] = j-1;
}

//calculate the total number of partitions
void updatePartitions(int nPart[], double partBoundX[], double partBoundY[], double L, double rad){
    nPart[0] = (int)(L/(2*rad));
    nPart[1] = (int)(L/(2*rad));
    for(int i=0; i<nPart[0]; i++)
        partBoundX[i] = i*2*rad;
    for(int i=0; i<nPart[1]; i++)
        partBoundY[i] = i*2*rad;
}

void projectIntoNewArea(particle *p, double L, double LOld) {
    for(int j=0; j<2; j++)
        p->position[j] = p->position[j]/LOld*L;
}

double repulsiveSpringPotential(double x1, double y1, double x2, double y2, double r1, double r2) {
    double epsilon0 = 1;
    double position1[2] = {x1,y1};
    double position2[2] = {x2,y2};
    
    double r[2];
    
    vectorSubtract(position2, position1, r);
    double R = vectorNorm(r);
    vectorNormalize(r, r);
    
    double epsilon = epsilon0;
    double sigma = r1+r2;
    
    if(R < sigma)
        return epsilon/2*pow(1-R/sigma, 2);
    else
        return 0;
}
/*
void repulsiveSpringForce(particle *p1, particle *p2, double L, int nPart[]) {
    double x1 = p1->position[0];
    double y1 = p1->position[1];
    double x2;
    double y2;
    
    double r1 = p1->rad;
    double r2 = p2->rad;
    
    if(p2->coord[0]==p1->coord[0]||abs(p2->coord[0]-p1->coord[0])==1)
        x2 = p2->position[0];
    else if(p2->coord[0]-p1->coord[0]==nPart[0]-1)
        x2 = p2->position[0]-L;
    else if(p2->coord[0]-p1->coord[0]==-(nPart[0]-1))
        x2 = p2->position[0]+L;
    else
        return;
    
    if(p2->coord[1]==p1->coord[1]||abs(p2->coord[1]-p1->coord[1])==1)
        y2 = p2->position[1];
    else if(p2->coord[1]-p1->coord[1]==nPart[1]-1)
        y2 = p2->position[1]-L;
    else if(p2->coord[1]-p1->coord[1]==-(nPart[1]-1))
        y2 = p2->position[1]+L;
    else
        return;
    
    if(!overlapQ(p1, p2, nPart, L))
        return;
    
    double delta = 1e-8;
    p1->force[0] = p1->force[0]+(repulsiveSpringPotential(x1-delta,y1,x2,y2,r1,r2)-repulsiveSpringPotential(x1+delta,y1,x2,y2,r1,r2))/2/delta;
    p1->force[1] = p1->force[1]+(repulsiveSpringPotential(x1,y1-delta,x2,y2,r1,r2)-repulsiveSpringPotential(x1,y1+delta,x2,y2,r1,r2))/2/delta;
    p2->force[0] = p2->force[0]+(repulsiveSpringPotential(x1,y1,x2-delta,y2,r1,r2)-repulsiveSpringPotential(x1,y1,x2+delta,y2,r1,r2))/2/delta;
    p2->force[1] = p2->force[1]+(repulsiveSpringPotential(x1,y1,x2,y2-delta,r1,r2)-repulsiveSpringPotential(x1,y1,x2,y2+delta,r1,r2))/2/delta;
   
}
*/
void repulsiveSpringForce(particle *p1, particle *p2, double L, int nPart[]) {
    double x1 = p1->position[0];
    double y1 = p1->position[1];
    double x2;
    double y2;
    
    double r1 = p1->rad;
    double r2 = p2->rad;
    
    if(p2->coord[0]==p1->coord[0]||abs(p2->coord[0]-p1->coord[0])==1)
        x2 = p2->position[0];
    else if(p2->coord[0]-p1->coord[0]==nPart[0]-1)
        x2 = p2->position[0]-L;
    else if(p2->coord[0]-p1->coord[0]==-(nPart[0]-1))
        x2 = p2->position[0]+L;
    else
        return;
    
    if(p2->coord[1]==p1->coord[1]||abs(p2->coord[1]-p1->coord[1])==1)
        y2 = p2->position[1];
    else if(p2->coord[1]-p1->coord[1]==nPart[1]-1)
        y2 = p2->position[1]-L;
    else if(p2->coord[1]-p1->coord[1]==-(nPart[1]-1))
        y2 = p2->position[1]+L;
    else
        return;
    
    if(!overlapQ(p1, p2, nPart, L))
        return;
    
    double epsilon = 1;
    double sigma = r1+r2;
    double R = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    double magnitude = epsilon/sigma*(1-R/sigma);
    p1->force[0] += magnitude/R*(x1-x2);
    p1->force[1] += magnitude/R*(y1-y2);
    p2->force[0] += magnitude/R*(x2-x1);
    p2->force[1] += magnitude/R*(y2-y1);
    
}

void addRepulsiveSpringForce (particle p[], int np, double L, int nPart[]) {
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++)
            repulsiveSpringForce(&(p[i]), &(p[j]), L, nPart);
}

double totalRepulsiveSpringEnergy(particle p[], int np, double L, int nPart[]) {
    double energy = 0;
    for(int i=0; i<np; i++)
        for(int j=i+1; j<np; j++) {
            double x1 = p[i].position[0];
            double y1 = p[i].position[1];
            double x2;
            double y2;
            
            if(p[j].coord[0]==p[i].coord[0]||abs(p[j].coord[0]-p[i].coord[0])==1)
                x2 = p[j].position[0];
            else if(p[j].coord[0]-p[i].coord[0]==nPart[0]-1)
                x2 = p[j].position[0]-L;
            else if(p[j].coord[0]-p[i].coord[0]==-(nPart[0]-1))
                x2 = p[j].position[0]+L;
            else
                continue;
            
            if(p[j].coord[1]==p[i].coord[1]||abs(p[j].coord[1]-p[i].coord[1])==1)
                y2 = p[j].position[1];
            else if(p[j].coord[1]-p[i].coord[1]==nPart[1]-1)
                y2 = p[j].position[1]-L;
            else if(p[j].coord[1]-p[i].coord[1]==-(nPart[1]-1))
                y2 = p[j].position[1]+L;
            else
                continue;
            
            energy = energy+repulsiveSpringPotential(x1, y1, x2, y2, p[i].rad, p[j].rad);
        }
    return energy;
}

void conjugateGradDescentIntegrate(particle *p, double dt, double L) {
    
    for (int i=0; i<2; i++) {
        p->position[i] = p->position[i] + (p->conjGrad[i])*dt;
    }
    while(p->position[0]<0)
        p->position[0] = p->position[0]+L;
    while(p->position[0]>L)
        p->position[0] = p->position[0]-L;
    while(p->position[1]<0)
        p->position[1] = p->position[1]+L;
    while(p->position[1]>L)
        p->position[1] = p->position[1]-L;
    
}

void calculateConjugateForce(particle p[], int np, double L, int nPart[]){
    //addGayBernePotentialForce(p, a, b, Rad, np);
    addRepulsiveSpringForce(p, np, L, nPart);
    double gammaNum = 0;
    double gammaDenom = 0;
    for(int i=0; i<np; i++)
        for(int j=0; j<2; j++){
            //gammaNum += (p[i].force[j]-p[i].oldForce[j])*p[i].force[j];
            gammaNum += p[i].force[j]*p[i].force[j];
            gammaDenom += p[i].oldForce[j]*p[i].oldForce[j];
        }
    for(int i=0; i<np; i++)
        for(int j=0; j<2; j++){
            p[i].conjGrad[j] = p[i].force[j] + gammaNum/gammaDenom*p[i].conjGrad[j];
            p[i].oldForce[j] = p[i].force[j];
        }
}

double energyForTimestep(particle p[], int np, double dt, double L, double partBoundX[], double partBoundY[], int nPart[]){
    for (int i=0; i<np; i++) {
        for(int j=0; j<2; j++)
            p[i].oldPosition[j] = p[i].position[j];
        for(int j=0; j<2; j++) {
            p[i].oldCoord[j] = p[i].coord[j];
        }
    }
    
    for(int i=0; i<np; i++){
        conjugateGradDescentIntegrate(&p[i], dt, L);
        partitionOneParticle(&p[i], partBoundX, partBoundY, nPart);
    }

    double energy = totalRepulsiveSpringEnergy(p, np, L, nPart);
    for (int i=0; i<np; i++) {
        for(int j=0; j<2; j++)
            p[i].position[j] = p[i].oldPosition[j];
        for(int j=0; j<2; j++) {
            p[i].coord[j] = p[i].oldCoord[j];
        }
    }
    return energy;
}

void newtonMove(particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[]){
    for(int i=0; i<np; i++){
        for(int j=0; j<2; j++)
            p[i].position[j] = p[i].position[j] + p[i].force[j];
        while(p[i].position[0]<0)
            p[i].position[0] = p[i].position[0]+L;
        while(p[i].position[0]>L)
            p[i].position[0] = p[i].position[0]-L;
        while(p[i].position[1]<0)
            p[i].position[1] = p[i].position[1]+L;
        while(p[i].position[1]>L)
            p[i].position[1] = p[i].position[1]-L;
        partitionOneParticle(&p[i], partBoundX, partBoundY, nPart);
    }
}

//find out the timestep that minimize the energy along the fixed conjugate gradient direction.
double goldenSearch(particle p[], int np, double rad, double L, double partBoundX[], double partBoundY[], int nPart[]) {
    double phi = 1.618033988749895;//goldenRatio+1
    double w = 0.3819660112501051;//1-goldenRatio
    double tol = MACHINE_EPSILON;
    //double currentEnergy = totalGayBernePotentialEnergy(p, np, a, b, Rad);
    double currentEnergy = totalRepulsiveSpringEnergy(p, np, L, nPart);
    
    double maxForce = 0;
    for(int i=0; i<np; i++)
        if(vectorNorm(p[i].conjGrad)>maxForce)
            maxForce = vectorNorm(p[i].conjGrad);
    double bracketLimit;
    double initialBracketGuess;
    if(maxForce>0) {
        bracketLimit = rad/maxForce;
        initialBracketGuess = 0.1*rad/maxForce;
    }
    else
        return 1;// return any finite real number
    
    double dtA, dtB, enA, enB;
    dtA = initialBracketGuess;
    enA = energyForTimestep(p, np, dtA, L, partBoundX, partBoundY, nPart);
    dtB = dtA*phi;
    enB = energyForTimestep(p, np, dtB, L, partBoundX, partBoundY, nPart);
    if(enA == 0)
        return dtA;
    do {
        if(enA >= currentEnergy) {
            if((enA-currentEnergy)/currentEnergy>=1e-16){
                dtB = dtA;
                enB = enA;
                dtA = dtA/phi;
                enA = energyForTimestep(p, np, dtA, L, partBoundX, partBoundY, nPart);
                if(enA == 0)
                    return dtA;
            } else
                return 0;
        } else if(enB <= enA) {
            dtA = dtB;
            enA = enB;
            dtB = dtB*phi;
            enB = energyForTimestep(p, np, dtB, L, partBoundX, partBoundY, nPart);
            if(enA == 0)
                return dtA;
            if(dtA > bracketLimit)
                return bracketLimit;
        }
    } while(enA>=currentEnergy || enB<=enA);
    //Now the initial bracket is (0, dtA, dtB). Find the minimum next.
    double dtx, dty;
    double enx, eny;
    dtx = dtA;
    enx = enA;
    dtA = 0;
    enA = currentEnergy;
    //(dtA, dtx, dty, dtB)
    do {
        if(dtB-dtx > dtx-dtA) {
            dty = dtx + w*(dtB-dtx);
            eny = energyForTimestep(p, np, dty, L, partBoundX, partBoundY, nPart);
        }
        else {
            dty = dtx;
            eny = enx;
            dtx = dtx - w*(dtx-dtA);
            enx = energyForTimestep(p, np, dtx, L, partBoundX, partBoundY, nPart);
        }
        if(enx<eny) {
            dtB = dty;
            enB = eny;
        } else {
            dtA = dtx;
            enA = enx;
            dtx = dty;
            enx = eny;
        }
    } while(dtB-dtA > tol*(dtx+dty));
    //printf("dt: %lf\n", dtx);
    return dtx;
}

void oneConjugateGradientDescentMoveGoldenSearch(particle p[], int np, double dt, double L, double partBoundX[], double partBoundY[], int nPart[]) {
    for(int i=0; i<np; i++){
        conjugateGradDescentIntegrate(&p[i], dt, L);
        partitionOneParticle(&p[i], partBoundX, partBoundY, nPart);
    }
}

void conjugateGradientDescentGoldenSearch(particle p[], int np, double rad, double L, double partBoundX[], double partBoundY[], int nPart[]) {
    double energy = totalRepulsiveSpringEnergy(p, np, L, nPart);
    double energyOld = energy;
    double dt;
    //addGayBernePotentialForce(p, a, b, Rad, np);
    addRepulsiveSpringForce(p, np, L, nPart);
    for(int i=0; i<np; i++) {
        for(int j=0; j<2; j++){
            p[i].conjGrad[j] = p[i].force[j];
            p[i].oldForce[j] = p[i].force[j];
        }
//        printf("%d: %.15lf, %.15lf\n", i, p[i].force[0], p[i].force[1]);
    }
    int conjGradCount = 0;
    while(1) {
        dt = goldenSearch(p, np, rad, L, partBoundX, partBoundY, nPart);
        if(dt == 0){//can't find minimum along conjugate direction
            if(conjGradCount==0){//can't find minimum even after reset conjugate direction, check directions from newton move
/*                int info;
                solveHessianSphere(p, np, L, nPart, partBoundX, partBoundY, &info);
                if(info==0) {
                    for(int j=0; j<np; j++)
                        for(int k=0; k<2; k++)
                            p[j].conjGrad[k] = p[j].force[k];
                    dt = goldenSearch(p, np, rad, L, partBoundX, partBoundY, nPart);
                    if(dt == 0) //cannot find minimum along the directions from newton method
                        break;
                    printf("Move along directions from Newton method!\n");
                    oneConjugateGradientDescentMoveGoldenSearch(p, np, dt, L, partBoundX, partBoundY, nPart);
                    energy = totalRepulsiveSpringEnergy(p, np, L, nPart);
                    //printf("%.25lf\n",energy);
                    if(energy/np<1e-16 || fabs(energy-energyOld)/energy<1e-16)
                        break;
                    energyOld = energy;
                    for(int j=0; j<np; j++)
                        for(int k=0; k<2; k++)
                            p[j].force[k] = 0;
                    addRepulsiveSpringForce(p, np, L, nPart);
                    for(int j=0; j<np; j++)
                        for(int k=0; k<2; k++){
                            p[j].conjGrad[k] = p[j].force[k];
                            p[j].oldForce[k] = p[j].force[k];
                        }
                    continue;
                }*/
                break;
            }
            else { //reset conjugate direction and try again
                conjGradCount = 0;
                for(int j=0; j<np; j++)
                    for(int k=0; k<2; k++){
                        p[j].conjGrad[k] = p[j].force[k];
                        p[j].oldForce[k] = p[j].force[k];
                    }
            }
        } else {
            /*
             FILE *conj = fopen("conjProj.dat","w");
             for(int i=0; i<np; i++) {
             fprintf(conj, "%.20lf %.20lf %.20lf %.20lf\n", p[i].conjGrad[0], p[i].conjGrad[1], p[i].conjGrad[2], p[i].conjGrad[3]);
             }
             fclose(conj);
             */
            oneConjugateGradientDescentMoveGoldenSearch(p, np, dt, L, partBoundX, partBoundY, nPart);
            conjGradCount++;
            if(conjGradCount>=100){
                int info;
                solveHessianSphere(p, np, L, nPart, partBoundX, partBoundY, &info);
                if(info==0) {
                    printf("Do one newton move!\n");
                    //newtonMove(p, np, a, b, L, nPart, partBoundX, partBoundY);
                    for(int j=0; j<np; j++)
                        for(int k=0; k<2; k++)
                            p[j].conjGrad[k] = p[j].force[k];
                    dt = goldenSearch(p, np, rad, L, partBoundX, partBoundY, nPart);
                    oneConjugateGradientDescentMoveGoldenSearch(p, np, dt, L, partBoundX, partBoundY, nPart);
                }
                else
                    printf("newton method failed!\n");
                conjGradCount = 0;
                energyOld = totalRepulsiveSpringEnergy(p, np, L, nPart);
                for(int j=0; j<np; j++)
                    for(int k=0; k<2; k++)
                        p[j].force[k] = 0;
                addRepulsiveSpringForce(p, np, L, nPart);
                for(int j=0; j<np; j++)
                    for(int k=0; k<2; k++){
                        p[j].conjGrad[k] = p[j].force[k];
                        p[j].oldForce[k] = p[j].force[k];
                    }
                continue;
            }
            
            energy = totalRepulsiveSpringEnergy(p, np, L, nPart);
            //printf("%d: %.25lf\n", conjGradCount, energy/np);
            if(energy/np<1e-16 || fabs(energy-energyOld)/energy<1e-16)
                break;
            energyOld = energy;
            for(int j=0; j<np; j++)
                for(int k=0; k<2; k++)
                    p[j].force[k] = 0;
            calculateConjugateForce(p, np, L, nPart);
        }
    }
    
}

void conjugateGradientDescentGoldenSearchSteps(particle p[], int np, double rad, double L, double partBoundX[], double partBoundY[], int nPart[], int steps) {
    double energy = totalRepulsiveSpringEnergy(p, np, L, nPart);
    double energyOld = energy;
    double dt;
    
    addRepulsiveSpringForce(p, np, L, nPart);
    for(int i=0; i<np; i++) {
        for(int j=0; j<2; j++){
            p[i].conjGrad[j] = p[i].force[j];
            p[i].oldForce[j] = p[i].force[j];
        }
        //        printf("%d: %.15lf, %.15lf\n", i, p[i].force[0], p[i].force[1]);
    }
    int conjGradCount = 0;
    int totalCount = 0;
    while(1) {
        dt = goldenSearch(p, np, rad, L, partBoundX, partBoundY, nPart);
        if(dt == 0){//can't find minimum along conjugate direction
            if(conjGradCount==0){
                break;
            }
            else { //reset conjugate direction and try again
                conjGradCount = 0;
                for(int j=0; j<np; j++)
                    for(int k=0; k<2; k++){
                        p[j].conjGrad[k] = p[j].force[k];
                        p[j].oldForce[k] = p[j].force[k];
                    }
            }
        } else {
            oneConjugateGradientDescentMoveGoldenSearch(p, np, dt, L, partBoundX, partBoundY, nPart);
            conjGradCount++;
            totalCount++;
            if(conjGradCount>=1000){
/*                int info;
                solveHessianSphere(p, np, L, nPart, partBoundX, partBoundY, &info);
                if(info==0) {
                    printf("Do one newton move!\n");
                    //newtonMove(p, np, a, b, L, nPart, partBoundX, partBoundY);
                    for(int j=0; j<np; j++)
                        for(int k=0; k<2; k++)
                            p[j].conjGrad[k] = p[j].force[k];
                    dt = goldenSearch(p, np, rad, L, partBoundX, partBoundY, nPart);
                    oneConjugateGradientDescentMoveGoldenSearch(p, np, dt, L, partBoundX, partBoundY, nPart);
                }
                else
                    printf("newton method failed!\n");
 */
                conjGradCount = 0;
                energyOld = totalRepulsiveSpringEnergy(p, np, L, nPart);
                for(int j=0; j<np; j++)
                    for(int k=0; k<2; k++)
                        p[j].force[k] = 0;
                addRepulsiveSpringForce(p, np, L, nPart);
                for(int j=0; j<np; j++)
                    for(int k=0; k<2; k++){
                        p[j].conjGrad[k] = p[j].force[k];
                        p[j].oldForce[k] = p[j].force[k];
                    }
                continue;
            }
            
            energy = totalRepulsiveSpringEnergy(p, np, L, nPart);
            //printf("%d: %.25lf\n", conjGradCount, energy/np);
            if(energy/np<1e-16 || fabs(energy-energyOld)/energyOld<1e-8 || totalCount>steps)
                break;
            energyOld = energy;
            for(int j=0; j<np; j++)
                for(int k=0; k<2; k++)
                    p[j].force[k] = 0;
            calculateConjugateForce(p, np, L, nPart);
        }
    }
    
}
