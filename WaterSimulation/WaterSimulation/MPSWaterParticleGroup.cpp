#include"MPSWaterParticleGroup.h"

void initializeParticlePositionAndVelocity_for2dim(void);
void initializeParticlePositionAndVelocity_for3dim(void);
void calculateConstantParameter(void);
void calculateNZeroAndLambda(void);
double weight(double distance, double re);
void calculateGravity(void);
void calculateViscosity(void);
void moveParticle(void);
void collision(void);
void calculatePressure(void);
void calculateParticleNumberDensity(void);
void setBoundaryCondition(void);
void setSourceTerm(void);
void setMatrix(void);
void exceptionalProcessingForBoundaryCondition(void);
void checkBoundaryCondition(void);
void increaseDiagonalTerm(void);
void solveSimultaniousEquationsByGaussEliminationMethod(void);
void removeNegativePressure(void);
void setMinimumPressure(void);
void calculatePressureGradient(void);
void moveParticleUsingPressureGradient(void);
static double AccelerationX[ARRAY_SIZE];
static double AccelerationY[ARRAY_SIZE];
static double AccelerationZ[ARRAY_SIZE];
static int    ParticleType[ARRAY_SIZE];
static float PositionX[ARRAY_SIZE];
static float PositionY[ARRAY_SIZE];
static float PositionZ[ARRAY_SIZE];
static double VelocityX[ARRAY_SIZE];
static double VelocityY[ARRAY_SIZE];
static double VelocityZ[ARRAY_SIZE];
static double Pressure[ARRAY_SIZE];
static double LastPressure[ARRAY_SIZE];
static double ParticleNumberDensity[ARRAY_SIZE];
static int    BoundaryCondition[ARRAY_SIZE];
static double SourceTerm[ARRAY_SIZE];
static int    FlagForCheckingBoundaryCondition[ARRAY_SIZE];
static double CoefficientMatrix[ARRAY_SIZE * ARRAY_SIZE];
static double MinimumPressure[ARRAY_SIZE];
int    FileNumber;
double Time;
int    NumberOfParticles;
double Re_forParticleNumberDensity, Re2_forParticleNumberDensity;
double Re_forGradient, Re2_forGradient;
double Re_forLaplacian, Re2_forLaplacian;
double N0_forParticleNumberDensity;
double N0_forGradient;
double N0_forLaplacian;
double Lambda;
double collisionDistance, collisionDistance2;
double FluidDensity;
int N0Count;


void initializeParticlePositionAndVelocity_for2dim(void) {

    int iX, iY;
    int nX, nY;
    double x, y, z;
    int i = 0;
    int flagOfParticleGeneration;

    nX = (int)(1.0 / PARTICLE_DISTANCE) + 5;
    nY = (int)(0.6 / PARTICLE_DISTANCE) + 5;
    for (iX = -4; iX < nX; iX++) {
        for (iY = -4; iY < nY; iY++) {
            x = PARTICLE_DISTANCE * (double)(iX);
            y = PARTICLE_DISTANCE * (double)(iY);
            z = 0.0;
            flagOfParticleGeneration = OFF;

            if (((x > -4.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS))) {  /* dummy wall region */
                ParticleType[i] = DUMMY_WALL;
                flagOfParticleGeneration = ON;
            }

            if (((x > -2.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 2.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.0 - 2.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS))) { /* wall region */
                ParticleType[i] = WALL;
                flagOfParticleGeneration = ON;
            }

            if (((x > -4.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.6 - 2.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS))) {  /* wall region */
                ParticleType[i] = WALL;
                flagOfParticleGeneration = ON;
            }

            if (((x > 0.0 + EPS) && (x <= 1.00 + EPS)) && (y > 0.0 + EPS)) {  /* empty region */
                flagOfParticleGeneration = OFF;
            }

            if (((x > 0.0 + EPS) && (x <= 0.25 + EPS)) && ((y > 0.0 + EPS) && (y <= 0.30 + EPS))) {  /* fluid region */
                ParticleType[i] = FLUID;
                flagOfParticleGeneration = ON;
            }

            if (flagOfParticleGeneration == ON) {
                PositionX[i] = x; PositionY[i] = y; PositionZ[i] = z;
                i++;
            }
        }
    }
    NumberOfParticles = i;
    for (i = 0; i < NumberOfParticles; i++) { VelocityX[i] = 0.0; VelocityY[i] = 0.0; VelocityZ[i] = 0.0; }
}

void initializeParticlePositionAndVelocity_for3dim(void) {
    int iX, iY, iZ;
    int nX, nY, nZ;
    double x, y, z;
    int i = 0;
    int flagOfParticleGeneration;

    nX = (int)(1.0 / PARTICLE_DISTANCE) + 5;
    nY = (int)(0.6 / PARTICLE_DISTANCE) + 5;
    nZ = (int)(0.3 / PARTICLE_DISTANCE) + 5;
    for (iX = -4; iX < nX; iX++) {
        for (iY = -4; iY < nY; iY++) {
            for (iZ = -4; iZ < nZ; iZ++) {
                x = PARTICLE_DISTANCE * iX;
                y = PARTICLE_DISTANCE * iY;
                z = PARTICLE_DISTANCE * iZ;
                flagOfParticleGeneration = OFF;

                /* dummy wall region */
                if ((((x > -4.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS))) && ((z > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) && (z <= 0.3 + 4.0 * PARTICLE_DISTANCE + EPS))) {
                    ParticleType[i] = DUMMY_WALL;
                    flagOfParticleGeneration = ON;
                }

                /* wall region */
                if ((((x > -2.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 2.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.0 - 2.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS))) && ((z > 0.0 - 2.0 * PARTICLE_DISTANCE + EPS) && (z <= 0.3 + 2.0 * PARTICLE_DISTANCE + EPS))) {
                    ParticleType[i] = WALL;
                    flagOfParticleGeneration = ON;
                }

                /* wall region */
                if ((((x > -4.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.6 - 2.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS))) && ((z > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) && (z <= 0.3 + 4.0 * PARTICLE_DISTANCE + EPS))) {
                    ParticleType[i] = WALL;
                    flagOfParticleGeneration = ON;
                }

                /* empty region */
                if ((((x > 0.0 + EPS) && (x <= 1.00 + EPS)) && (y > 0.0 + EPS)) && ((z > 0.0 + EPS) && (z <= 0.3 + EPS))) {
                    flagOfParticleGeneration = OFF;
                }

                /* fluid region */
                if ((((x > 0.0 + EPS) && (x <= 0.25 + EPS)) && ((y > 0.0 + EPS) && (y < 0.5 + EPS))) && ((z > 0.0 + EPS) && (z <= 0.3 + EPS))) {
                    ParticleType[i] = FLUID;
                    flagOfParticleGeneration = ON;
                }

                if (flagOfParticleGeneration == ON) {
                    PositionX[i] = x;
                    PositionY[i] = y;
                    PositionZ[i] = z;
                    i++;
                }
            }
        }
    }
    NumberOfParticles = i;
    for (i = 0; i < NumberOfParticles; i++) { VelocityX[i] = 0.0; VelocityY[i] = 0.0; VelocityZ[i] = 0.0; }
}

void calculateConstantParameter(void) {

    Re_forParticleNumberDensity = RADIUS_FOR_NUMBER_DENSITY;
    Re_forGradient = RADIUS_FOR_GRADIENT;
    Re_forLaplacian = RADIUS_FOR_LAPLACIAN;
    Re2_forParticleNumberDensity = Re_forParticleNumberDensity * Re_forParticleNumberDensity;
    Re2_forGradient = Re_forGradient * Re_forGradient;
    Re2_forLaplacian = Re_forLaplacian * Re_forLaplacian;
    calculateNZeroAndLambda();
    FluidDensity = FLUID_DENSITY;
    collisionDistance = COLLISION_DISTANCE;
    collisionDistance2 = collisionDistance * collisionDistance;
    FileNumber = 0;
    Time = 0.0;
}

void calculateNZeroAndLambda(void) {
    int iX, iY, iZ;
    int iZ_start, iZ_end;
    double xj, yj, zj, distance, distance2;
    double xi, yi, zi;

    if (DIM == 2) {
        iZ_start = 0; iZ_end = 1;
    }
    else {
        iZ_start = -4; iZ_end = 5;
    }

    N0_forParticleNumberDensity = 0.0;
    N0_forGradient = 0.0;
    N0_forLaplacian = 0.0;
    Lambda = 0.0;
    xi = 0.0;  yi = 0.0;  zi = 0.0;

    N0Count = 0;

    for (iX = -4; iX < 5; iX++) {
        for (iY = -4; iY < 5; iY++) {
            for (iZ = iZ_start; iZ < iZ_end; iZ++) {
                if (((iX == 0) && (iY == 0)) && (iZ == 0))continue;
                xj = PARTICLE_DISTANCE * (double)(iX);
                yj = PARTICLE_DISTANCE * (double)(iY);
                zj = PARTICLE_DISTANCE * (double)(iZ);
                distance2 = (xj - xi) * (xj - xi) + (yj - yi) * (yj - yi) + (zj - zi) * (zj - zi);
                distance = sqrt(distance2);
                if (weight(distance, Re_forParticleNumberDensity) != 0)
                {
                    N0Count++;
                }
                N0_forParticleNumberDensity += weight(distance, Re_forParticleNumberDensity);
                N0_forGradient += weight(distance, Re_forGradient);
                N0_forLaplacian += weight(distance, Re_forLaplacian);
                Lambda += distance2 * weight(distance, Re_forLaplacian);
            }
        }
    }
    Lambda = Lambda / N0_forLaplacian;
}


double weight(double distance, double re) {
    double weightIJ;

    if (distance >= re) {
        weightIJ = 0.0;
    }
    else {
        weightIJ = (re / distance) - 1.0;
    }
    return weightIJ;
}

void calculateGravity(void) {
    int i;
#pragma omp parallel for
    for (i = 0; i < NumberOfParticles; i++) {
        if (ParticleType[i] == FLUID) {
            AccelerationX[i] = GRAVITY_X;
            AccelerationY[i] = GRAVITY_Y;
            AccelerationZ[i] = GRAVITY_Z;
        }
        else {
            AccelerationX[i] = 0.0;
            AccelerationY[i] = 0.0;
            AccelerationZ[i] = 0.0;
        }
    }
}

void calculateViscosity(void) {
    double a;

    a = (KINEMATIC_VISCOSITY) * (2.0 * DIM) / (N0_forLaplacian * Lambda);

#pragma omp parallel for
    for (int i = 0; i < NumberOfParticles; i++) {
        if (ParticleType[i] != FLUID) continue;
        double viscosityTermX, viscosityTermY, viscosityTermZ;
        viscosityTermX = 0.0;  viscosityTermY = 0.0;  viscosityTermZ = 0.0;
        //#pragma omp parallel for reduction(+:viscosityTermX) reduction(+:viscosityTermY) reduction(+:viscosityTermZ)
        for (int j = 0; j < NumberOfParticles; j++) {
            if ((j == i) || (ParticleType[j] == GHOST)) continue;
            double xij, yij, zij;
            double distance, distance2;
            xij = PositionX[j] - PositionX[i];
            yij = PositionY[j] - PositionY[i];
            zij = PositionZ[j] - PositionZ[i];
            distance2 = (xij * xij) + (yij * yij) + (zij * zij);
            distance = sqrt(distance2);
            if (distance < Re_forLaplacian) {
                double w;
                w = weight(distance, Re_forLaplacian);
                viscosityTermX += (VelocityX[j] - VelocityX[i]) * w;
                viscosityTermY += (VelocityY[j] - VelocityY[i]) * w;
                viscosityTermZ += (VelocityZ[j] - VelocityZ[i]) * w;
            }
        }
        viscosityTermX = viscosityTermX * a;
        viscosityTermY = viscosityTermY * a;
        viscosityTermZ = viscosityTermZ * a;
        AccelerationX[i] += viscosityTermX;
        AccelerationY[i] += viscosityTermY;
        AccelerationZ[i] += viscosityTermZ;
    }
}

void moveParticle(void) {
#pragma omp parallel for
    for (int i = 0; i < NumberOfParticles; i++) {
        if (ParticleType[i] == FLUID) {
            VelocityX[i] += AccelerationX[i] * DT;
            VelocityY[i] += AccelerationY[i] * DT;
            VelocityZ[i] += AccelerationZ[i] * DT;

            PositionX[i] += VelocityX[i] * DT;
            PositionY[i] += VelocityY[i] * DT;
            PositionZ[i] += VelocityZ[i] * DT;
        }
        AccelerationX[i] = 0.0;
        AccelerationY[i] = 0.0;
        AccelerationZ[i] = 0.0;
    }
}

void collision(void) {
    int    i, j;
    double xij, yij, zij;
    double distance, distance2;
    double forceDT; /* forceDT is the impulse of collision between particles */
    double mi, mj;
    double velocity_ix, velocity_iy, velocity_iz;
    double e = COEFFICIENT_OF_RESTITUTION;
    static double VelocityAfterCollisionX[ARRAY_SIZE];
    static double VelocityAfterCollisionY[ARRAY_SIZE];
    static double VelocityAfterCollisionZ[ARRAY_SIZE];

    for (i = 0; i < NumberOfParticles; i++) {
        VelocityAfterCollisionX[i] = VelocityX[i];
        VelocityAfterCollisionY[i] = VelocityY[i];
        VelocityAfterCollisionZ[i] = VelocityZ[i];
    }

//#pragma omp parallel for
    for (i = 0; i < NumberOfParticles; i++) {
        if (ParticleType[i] == FLUID) {
            mi = FluidDensity;
            velocity_ix = VelocityX[i];
            velocity_iy = VelocityY[i];
            velocity_iz = VelocityZ[i];
            for (j = 0; j < NumberOfParticles; j++) {
                if ((j == i) || (ParticleType[j] == GHOST)) continue;
                xij = PositionX[j] - PositionX[i];
                yij = PositionY[j] - PositionY[i];
                zij = PositionZ[j] - PositionZ[i];
                distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                if (distance2 < collisionDistance2) {
                    distance = sqrt(distance2);
                    forceDT = (velocity_ix - VelocityX[j]) * (xij / distance)
                        + (velocity_iy - VelocityY[j]) * (yij / distance)
                        + (velocity_iz - VelocityZ[j]) * (zij / distance);
                    if (forceDT > 0.0) {
                        mj = FluidDensity;
                        forceDT *= (1.0 + e) * mi * mj / (mi + mj);
                        velocity_ix -= (forceDT / mi) * (xij / distance);
                        velocity_iy -= (forceDT / mi) * (yij / distance);
                        velocity_iz -= (forceDT / mi) * (zij / distance);
                        /*
                        if(j>i){ fprintf(stderr,"WARNING: Collision occured between %d and %d particles.\n",i,j); }
                        */
                    }
                }
            }
            VelocityAfterCollisionX[i] = velocity_ix;
            VelocityAfterCollisionY[i] = velocity_iy;
            VelocityAfterCollisionZ[i] = velocity_iz;
        }
    }
    for (i = 0; i < NumberOfParticles; i++) {
        if (ParticleType[i] == FLUID) {
            PositionX[i] += (VelocityAfterCollisionX[i] - VelocityX[i]) * DT;
            PositionY[i] += (VelocityAfterCollisionY[i] - VelocityY[i]) * DT;
            PositionZ[i] += (VelocityAfterCollisionZ[i] - VelocityZ[i]) * DT;
            VelocityX[i] = VelocityAfterCollisionX[i];
            VelocityY[i] = VelocityAfterCollisionY[i];
            VelocityZ[i] = VelocityAfterCollisionZ[i];
        }
    }
}

void calculatePressure(void) {
    calculateParticleNumberDensity();
    setBoundaryCondition();
    setSourceTerm();
    setMatrix();
    solveSimultaniousEquationsByGaussEliminationMethod();
    removeNegativePressure();
    setMinimumPressure();
}

void calculateParticleNumberDensity(void) {
    //int    i, j;

#pragma omp parallel for
    for (int i = 0; i < NumberOfParticles; i++) {
        ParticleNumberDensity[i] = 0.0;
        if (ParticleType[i] == GHOST) continue;
        for (int j = 0; j < NumberOfParticles; j++) {
            if ((j == i) || (ParticleType[j] == GHOST)) continue;
            double w;
            double xij, yij, zij;
            double distance, distance2;
            xij = PositionX[j] - PositionX[i];
            yij = PositionY[j] - PositionY[i];
            zij = PositionZ[j] - PositionZ[i];
            distance2 = (xij * xij) + (yij * yij) + (zij * zij);
            distance = sqrt(distance2);
            w = weight(distance, Re_forParticleNumberDensity);
            ParticleNumberDensity[i] += w;
        }
    }
}

void setBoundaryCondition(void) {
    
    double n0 = N0_forParticleNumberDensity;
    double beta = THRESHOLD_RATIO_OF_NUMBER_DENSITY;

#pragma omp parallel for
    for (int i = 0; i < NumberOfParticles; i++) {
        //int temp = ParticleType[i];
        float count = 0;
        for (int j = 0; j < NumberOfParticles; j++)
        {
            double distance, distance2;
            double xij, yij, zij;
            xij = PositionX[j] - PositionX[i];
            yij = PositionY[j] - PositionY[i];
            zij = PositionZ[j] - PositionZ[i];
            distance2 = (xij * xij) + (yij * yij) + (zij * zij);
            distance = sqrt(distance2);
            if (distance < Re_forParticleNumberDensity)
                count++;
        }


        if (ParticleType[i] == GHOST || ParticleType[i] == DUMMY_WALL) {
            BoundaryCondition[i] = GHOST_OR_DUMMY;
        }
        else if (count < 0.95 * N0Count) {
            BoundaryCondition[i] = SURFACE_PARTICLE;
            //ParticleType[i] = SURFACE;
            //count < 0.95 * N0Count
            //ParticleNumberDensity[i] < beta * n0
        }
        else {
            //ParticleType[i] = temp;
            BoundaryCondition[i] = INNER_PARTICLE;
        }
    }
}

void setSourceTerm(void) {
    double n0 = N0_forParticleNumberDensity;
    double gamma = RELAXATION_COEFFICIENT_FOR_PRESSURE;

#pragma omp parallel for
    for (int i = 0; i < NumberOfParticles; i++) {
        SourceTerm[i] = 0.0;
        if (ParticleType[i] == GHOST || ParticleType[i] == DUMMY_WALL) continue;
        if (BoundaryCondition[i] == INNER_PARTICLE) {
            SourceTerm[i] = gamma * (1.0 / (DT * DT)) * ((ParticleNumberDensity[i] - n0) / n0);
        }
        else if (BoundaryCondition[i] == SURFACE_PARTICLE) {
            SourceTerm[i] = 0.0;
        }
    }
}

void setMatrix(void) {



    double n0 = N0_forLaplacian;
    //int    i, j;
    double a;
    int n = NumberOfParticles;

    for (int i = 0; i < NumberOfParticles; i++) {
        for (int j = 0; j < NumberOfParticles; j++) {
            CoefficientMatrix[i * n + j] = 0.0;
        }
    }

    a = 2.0 * DIM / (n0 * Lambda);

#pragma omp parallel for
    for (int i = 0; i < NumberOfParticles; i++) {
        if (BoundaryCondition[i] != INNER_PARTICLE) continue;
        for (int j = 0; j < NumberOfParticles; j++) {
            if ((j == i) || (BoundaryCondition[j] == GHOST_OR_DUMMY)) continue;
            double coefficientIJ;
            double distance, distance2;
            double xij, yij, zij;
            xij = PositionX[j] - PositionX[i];
            yij = PositionY[j] - PositionY[i];
            zij = PositionZ[j] - PositionZ[i];
            distance2 = (xij * xij) + (yij * yij) + (zij * zij);
            distance = sqrt(distance2);
            if (distance >= Re_forLaplacian)continue;
            coefficientIJ = a * weight(distance, Re_forLaplacian) / FluidDensity;
            CoefficientMatrix[i * n + j] = (-1.0) * coefficientIJ;
            CoefficientMatrix[i * n + i] += coefficientIJ;
        }
        CoefficientMatrix[i * n + i] += (COMPRESSIBILITY) / (DT * DT);
    }
    //exceptionalProcessingForBoundaryCondition();
}

void exceptionalProcessingForBoundaryCondition(void) {
    /* If tere is no Dirichlet boundary condition on the fluid,
       increase the diagonal terms of the matrix for an exception. This allows us to solve the matrix without Dirichlet boundary conditions. */
    checkBoundaryCondition();
    increaseDiagonalTerm();
}

void checkBoundaryCondition(void) {
    int i, j, count;
    double xij, yij, zij, distance2;

    for (i = 0; i < NumberOfParticles; i++) {
        if (BoundaryCondition[i] == GHOST_OR_DUMMY) {
            FlagForCheckingBoundaryCondition[i] = GHOST_OR_DUMMY;
        }
        else if (BoundaryCondition[i] == SURFACE_PARTICLE) {
            FlagForCheckingBoundaryCondition[i] = DIRICHLET_BOUNDARY_IS_CONNECTED;
        }
        else {
            FlagForCheckingBoundaryCondition[i] = DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
        }
    }

    do {
        count = 0;
        for (i = 0; i < NumberOfParticles; i++) {
            if (FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_CONNECTED) {
                for (j = 0; j < NumberOfParticles; j++) {
                    if (j == i) continue;
                    if ((ParticleType[j] == GHOST) || (ParticleType[j] == DUMMY_WALL)) continue;
                    if (FlagForCheckingBoundaryCondition[j] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) {
                        xij = PositionX[j] - PositionX[i];
                        yij = PositionY[j] - PositionY[i];
                        zij = PositionZ[j] - PositionZ[i];
                        distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                        if (distance2 >= Re2_forLaplacian)continue;
                        FlagForCheckingBoundaryCondition[j] = DIRICHLET_BOUNDARY_IS_CONNECTED;
                    }
                }
                FlagForCheckingBoundaryCondition[i] = DIRICHLET_BOUNDARY_IS_CHECKED;
                count++;
            }
        }
    } while (count != 0); /* This procedure is repeated until the all fluid or wall particles (which have Dirhchlet boundary condition in the particle group) are in the state of "DIRICHLET_BOUNDARY_IS_CHECKED".*/

    for (i = 0; i < NumberOfParticles; i++) {
        if (FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) {
            fprintf(stderr, "WARNING: There is no dirichlet boundary condition for %d-th particle.\n", i);
        }
    }
}

void increaseDiagonalTerm(void) {
    int i;
    int n = NumberOfParticles;

    for (i = 0; i < n; i++) {
        if (FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED) {
            CoefficientMatrix[i * n + i] = 2.0 * CoefficientMatrix[i * n + i];
        }
    }
}

void solveSimultaniousEquationsByGaussEliminationMethod(void) {
    //int    i, j, k;


    int    n = NumberOfParticles;

#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        Pressure[i] = 0.0;
    }

#pragma omp parallel for
    for (int i = 0; i < n - 1; i++) {
        if (BoundaryCondition[i] != INNER_PARTICLE) continue;
        for (int j = i + 1; j < n; j++) {
            if (BoundaryCondition[j] == GHOST_OR_DUMMY) continue;
            double c;
            c = CoefficientMatrix[j * n + i] / CoefficientMatrix[i * n + i];
            for (int k = i + 1; k < n; k++) {
                CoefficientMatrix[j * n + k] -= c * CoefficientMatrix[i * n + k];
            }
            SourceTerm[j] -= c * SourceTerm[i];
        }
    }

#pragma omp parallel for
    for (int i = n - 1; i >= 0; i--) {
        if (BoundaryCondition[i] != INNER_PARTICLE) continue;
        double sumOfTerms;
        sumOfTerms = 0.0;
        for (int j = i + 1; j < n; j++) {
            if (BoundaryCondition[j] == GHOST_OR_DUMMY) continue;
            sumOfTerms += CoefficientMatrix[i * n + j] * Pressure[j];
        }
        Pressure[i] = (SourceTerm[i] - sumOfTerms) / CoefficientMatrix[i * n + i];


    }

    //for (int i = 0; i < NumberOfParticles; i++)
    //{
    //    if (BoundaryCondition[i] == INNER_PARTICLE && ParticleType[i] == FLUID)
    //    {
    //        cout << Pressure[i] << endl;
    //    }
    //}

}

void removeNegativePressure(void) {
    int i;

    for (i = 0; i < NumberOfParticles; i++) {
        if (Pressure[i] < 0.0)Pressure[i] = 0.0;
    }
}

void setMinimumPressure(void) {
    double xij, yij, zij, distance2;
    int i, j;

    for (i = 0; i < NumberOfParticles; i++) {
        if (ParticleType[i] == GHOST || ParticleType[i] == DUMMY_WALL)continue;
        MinimumPressure[i] = Pressure[i];
        for (j = 0; j < NumberOfParticles; j++) {
            if ((j == i) || (ParticleType[j] == GHOST)) continue;
            if (ParticleType[j] == DUMMY_WALL) continue;
            xij = PositionX[j] - PositionX[i];
            yij = PositionY[j] - PositionY[i];
            zij = PositionZ[j] - PositionZ[i];
            distance2 = (xij * xij) + (yij * yij) + (zij * zij);
            if (distance2 >= Re2_forGradient)continue;
            if (MinimumPressure[i] > Pressure[j]) {
                MinimumPressure[i] = Pressure[j];
            }
        }
    }
}

void calculatePressureGradient(void) {
    int    i, j;
    double gradient_x, gradient_y, gradient_z;
    double xij, yij, zij;
    double distance, distance2;
    double w, pij;
    double a;

    a = DIM / N0_forGradient;
    for (i = 0; i < NumberOfParticles; i++) {
        if (ParticleType[i] != FLUID) continue;
        gradient_x = 0.0;  gradient_y = 0.0;  gradient_z = 0.0;
        for (j = 0; j < NumberOfParticles; j++) {
            if (j == i) continue;
            if (ParticleType[j] == GHOST) continue;
            if (ParticleType[j] == DUMMY_WALL) continue;
            xij = PositionX[j] - PositionX[i];
            yij = PositionY[j] - PositionY[i];
            zij = PositionZ[j] - PositionZ[i];
            distance2 = (xij * xij) + (yij * yij) + (zij * zij);
            distance = sqrt(distance2);
            if (distance < Re_forGradient) {
                w = weight(distance, Re_forGradient);
                pij = (Pressure[j] - MinimumPressure[i]) / distance2;
                gradient_x += xij * pij * w;
                gradient_y += yij * pij * w;
                gradient_z += zij * pij * w;
            }
        }
        gradient_x *= a;
        gradient_y *= a;
        gradient_z *= a;
        AccelerationX[i] = (-1.0) * gradient_x / FluidDensity;
        AccelerationY[i] = (-1.0) * gradient_y / FluidDensity;
        AccelerationZ[i] = (-1.0) * gradient_z / FluidDensity;
    }
}

void moveParticleUsingPressureGradient(void) {
    int i;

    for (i = 0; i < NumberOfParticles; i++) {
        if (ParticleType[i] == FLUID) {
            VelocityX[i] += AccelerationX[i] * DT;
            VelocityY[i] += AccelerationY[i] * DT;
            VelocityZ[i] += AccelerationZ[i] * DT;

            PositionX[i] += AccelerationX[i] * DT * DT;
            PositionY[i] += AccelerationY[i] * DT * DT;
            PositionZ[i] += AccelerationZ[i] * DT * DT;
        }
        AccelerationX[i] = 0.0;
        AccelerationY[i] = 0.0;
        AccelerationZ[i] = 0.0;
    }
}

void MPSWaterParticleGroup::InitParticles()
{
    printf("\n*** START MPS-SIMULATION ***\n");
    if (DIM == 2) {
        initializeParticlePositionAndVelocity_for2dim();
    }
    else {
        initializeParticlePositionAndVelocity_for3dim();
    }
    calculateConstantParameter();

    printf("*** END ***\n\n");

    //calculateGravity();
    //calculateViscosity();
    ////moveParticle();
    ////collision();
    //calculatePressure();
    //calculatePressureGradient();
    //moveParticleUsingPressureGradient();
}

void MPSWaterParticleGroup::Modeling()
{
	vector<vec3> posArray;
	for (int i = 0; i < particles.size(); i++)
	{
		posArray.push_back(particles[i].position);
	}

	vector<float> verticesInfo;
	bool provideNormal;
	bool provideTex;
	marchingCube = new MarchingCube();
	marchingCube->GetMeshData(posArray, verticesInfo, provideNormal, provideTex);

	shaderData->InitVertexBuffer(verticesInfo, provideNormal, provideTex);
}

void MPSWaterParticleGroup::Update(float dt)
{
    //cout <<  1 / dt << endl;
    //int iTimeStep = 0;
    shaderData->UpdateMatrix(transformation);
	//return;

    calculateGravity();
    calculateViscosity();
    moveParticle();
    collision();
    calculatePressure();
    calculatePressureGradient();
    moveParticleUsingPressureGradient();
    InitBufferData();
}

void MPSWaterParticleGroup::Draw()
{
	// GL_POINTS
	shaderData->SetDrawType(GL_POINTS);
	renderer->Render(shaderData);
}

void MPSWaterParticleGroup::InitBufferData()
{
	vector<float> data;

    for (int i = 0; i < NumberOfParticles; i++) 
    {
        vec3 color = vec3(0);
        switch (BoundaryCondition[i])
        {
        case INNER_PARTICLE:
            if (ParticleType[i] == FLUID)
                color = vec3(0, 0, 255);
            else
                //color = vec3(0);
                color = vec3(0, 255, 0);
            break;
        case SURFACE_PARTICLE:
            if (ParticleType[i] == FLUID)
                color = vec3(255, 0, 0);
            else
                color = vec3(0, 255, 0);
            break;
        case GHOST_OR_DUMMY:
            color = vec3(0);
            //color = vec3(0, 255, 255);
            break;
        //case WALL:
        //    color = vec3(0, 255, 0);
        //    break;
        default:
            break;
        }
        data.push_back(PositionX[i]);
        data.push_back(PositionY[i]);
        data.push_back(PositionZ[i]);
        data.push_back(color.x);
        data.push_back(color.y);
        data.push_back(color.z);
    }


	//for (int i = 0; i < particles.size(); i++)
	//{
	//	vec3 color = vec3(0);
	//	switch (particles[i].particleType)
	//	{
	//	case FLUID:
	//		color = vec3(0, 0, 255);
	//		break;
	//	case SURFACE:
	//		color = vec3(255, 0, 0);
	//		break;
	//	case DUMMY:
	//		color = vec3(0);
	//		//color = vec3(0, 255, 255);
	//		break;
	//	case WALL:
	//		color = vec3(0, 255, 0);
	//		break;
	//	default:
	//		break;
	//	}

	//	data.push_back(particles[i].position.x);
	//	data.push_back(particles[i].position.y);
	//	data.push_back(particles[i].position.z);
	//	data.push_back(color.x);
	//	data.push_back(color.y);
	//	data.push_back(color.z);
	//}

	shaderData->drawUnitNumber = data.size() / 6;
	dynamic_cast<VCShaderData*>(shaderData)->InitVertexBuffer(data, GL_DYNAMIC_DRAW);
}

void MPSWaterParticleGroup::UpdateBufferData()
{

}
