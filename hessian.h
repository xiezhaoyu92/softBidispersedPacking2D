//
//  hessian.h
//  softEllipsoid2D
//
//  Created by Zhaoyu Xie on 11/19/19.
//  Copyright © 2019 Zhaoyu Xie. All rights reserved.
//

#ifndef hessian_h
#define hessian_h

#include <stdio.h>
#include "particles.h"

void solveHessianSphere(particle p[], int np, double L, int nPart[], double partBoundX[], double partBoundY[], int *info);

#endif /* hessian_h */
