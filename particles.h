//
//  particles.h
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#ifndef particles_h
#define particles_h

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-15
#endif

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#include <stdio.h>

typedef struct {
    double position[2];
    double oldPosition[2];
    double postPreRelaxPosition[2];
    
    double force[2];
    double forceStoc[2];
    double oldForce[2];
    double conjGrad[2];
    double rad;
    int oldCoord[2];
    int coord[2];
    double postPreRelaxCoord[2];
} particle;

int overlapQ(particle *p1, particle *p2, int nPart[], double L);
void gradDescentConstrainedIntegrate(particle *p, double dt, double L);
void addStochasticForce(double D, double rad, particle *p);
void partitionOneParticle(particle *p, double partBoundX[], double partBoundY[], int nPart[]);
void updatePartitions(int nPart[], double partBoundX[], double partBoundY[], double L, double rad);
void projectIntoNewArea(particle *p, double L, double LOld);

double repulsiveSpringPotential(double x1, double y1, double x2, double y2, double r1, double r2);
void addRepulsiveSpringForce (particle p[], int np, double L, int nPart[]);
double totalRepulsiveSpringEnergy(particle p[], int np, double L, int nPart[]);

void conjugateGradientDescentGoldenSearch(particle p[], int np, double rad, double L, double partBoundX[], double partBoundY[], int nPart[]);
void conjugateGradientDescentGoldenSearchSteps(particle p[], int np, double rad, double L, double partBoundX[], double partBoundY[], int nPart[], int steps);


#endif /* particles_h */
