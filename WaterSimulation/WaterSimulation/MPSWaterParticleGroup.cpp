#define _CRT_SECURE_NO_WARNINGS

#include"MPSWaterParticleGroup.h"

/*=====================================================================
  mps.c

    Sample program of the MPS method

    Moving Particle Semi-implicit Method, Academic Press, 2018
    ISBN: 9780128127797

=======================================================================*/
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>



#define DIM                  2
#define PARTICLE_DISTANCE    0.025
#define DT                   0.001
#define OUTPUT_INTERVAL      20
//
///* for three-dimensional simulation */
///*
//#define DIM                  3
//#define PARTICLE_DISTANCE    0.075
//#define DT                   0.003
//#define OUTPUT_INTERVAL      2
//*/

//#define DIM                  3
//#define PARTICLE_DISTANCE    0.075
//#define DT                   0.003
//#define OUTPUT_INTERVAL      2

#define ARRAY_SIZE           5000
#define FINISH_TIME          2.0
#define KINEMATIC_VISCOSITY  (1.0E-6)
#define FLUID_DENSITY        1000.0 
#define GRAVITY_X  0.0      
#define GRAVITY_Y  -9.8
#define GRAVITY_Z  0.0      
#define RADIUS_FOR_NUMBER_DENSITY  (2.1*PARTICLE_DISTANCE) 
#define RADIUS_FOR_GRADIENT        (2.1*PARTICLE_DISTANCE) 
#define RADIUS_FOR_LAPLACIAN       (3.1*PARTICLE_DISTANCE) 
#define COLLISION_DISTANCE         (0.5*PARTICLE_DISTANCE)
#define THRESHOLD_RATIO_OF_NUMBER_DENSITY  0.97   
#define COEFFICIENT_OF_RESTITUTION 0.2
#define COMPRESSIBILITY (0.45E-9)
#define EPS             (0.01 * PARTICLE_DISTANCE)     
#define ON              1
#define OFF             0
#define RELAXATION_COEFFICIENT_FOR_PRESSURE 0.2
#define GHOST  -1
#define FLUID   0
#define SURFACE 1
#define WALL    2
#define DUMMY_WALL  3
#define GHOST_OR_DUMMY  -1
#define SURFACE_PARTICLE 1    
#define INNER_PARTICLE   0      
#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0 
#define DIRICHLET_BOUNDARY_IS_CONNECTED     1 
#define DIRICHLET_BOUNDARY_IS_CHECKED       2 
void initializeParticlePositionAndVelocity_for2dim(void);
void initializeParticlePositionAndVelocity_for3dim(void);
void calculateConstantParameter(void);
void calculateNZeroAndLambda(void);
double weight(double distance, double re);
void mainLoopOfSimulation(void);
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
void writeData_inProfFormat(void);
void writeData_inVtuFormat(void);
static double AccelerationX[ARRAY_SIZE];
static double AccelerationY[ARRAY_SIZE];
static double AccelerationZ[ARRAY_SIZE];
static int    ParticleType[ARRAY_SIZE];
static double PositionX[ARRAY_SIZE];
static double PositionY[ARRAY_SIZE];
static double PositionZ[ARRAY_SIZE];
static double VelocityX[ARRAY_SIZE];
static double VelocityY[ARRAY_SIZE];
static double VelocityZ[ARRAY_SIZE];
static double Pressure[ARRAY_SIZE];
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

//int main(void) {
//
//    printf("\n*** START MPS-SIMULATION ***\n");
//    if (DIM == 2) {
//        initializeParticlePositionAndVelocity_for2dim();
//    }
//    else {
//        initializeParticlePositionAndVelocity_for3dim();
//    }
//    calculateConstantParameter();
//    mainLoopOfSimulation();
//    printf("*** END ***\n\n");
//    return 0;
//
//}


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

            if (((x > 0.0 + EPS) && (x <= 0.25 + EPS)) && ((y > 0.0 + EPS) && (y <= 0.50 + EPS))) {  /* fluid region */
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

    for (iX = -4; iX < 5; iX++) {
        for (iY = -4; iY < 5; iY++) {
            for (iZ = iZ_start; iZ < iZ_end; iZ++) {
                if (((iX == 0) && (iY == 0)) && (iZ == 0))continue;
                xj = PARTICLE_DISTANCE * (double)(iX);
                yj = PARTICLE_DISTANCE * (double)(iY);
                zj = PARTICLE_DISTANCE * (double)(iZ);
                distance2 = (xj - xi) * (xj - xi) + (yj - yi) * (yj - yi) + (zj - zi) * (zj - zi);
                distance = sqrt(distance2);
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


void mainLoopOfSimulation(void) {
    int iTimeStep = 0;

    writeData_inVtuFormat();
    writeData_inProfFormat();

    while (1) {
        calculateGravity();
        calculateViscosity();
        moveParticle();
        collision();
        calculatePressure();
        calculatePressureGradient();
        moveParticleUsingPressureGradient();
        iTimeStep++;
        Time += DT;
        if ((iTimeStep % OUTPUT_INTERVAL) == 0) {
            printf("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticles: %d\n", iTimeStep, Time, NumberOfParticles);
            writeData_inVtuFormat();
            writeData_inProfFormat();
        }
        if (Time >= FINISH_TIME) { break; }
    }
}


void calculateGravity(void) {
    int i;

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
    int i, j;
    double viscosityTermX, viscosityTermY, viscosityTermZ;
    double distance, distance2;
    double w;
    double xij, yij, zij;
    double a;

    a = (KINEMATIC_VISCOSITY) * (2.0 * DIM) / (N0_forLaplacian * Lambda);
    for (i = 0; i < NumberOfParticles; i++) {
        if (ParticleType[i] != FLUID) continue;
        viscosityTermX = 0.0;  viscosityTermY = 0.0;  viscosityTermZ = 0.0;

        for (j = 0; j < NumberOfParticles; j++) {
            if ((j == i) || (ParticleType[j] == GHOST)) continue;
            xij = PositionX[j] - PositionX[i];
            yij = PositionY[j] - PositionY[i];
            zij = PositionZ[j] - PositionZ[i];
            distance2 = (xij * xij) + (yij * yij) + (zij * zij);
            distance = sqrt(distance2);
            if (distance < Re_forLaplacian) {
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
    int i;

    for (i = 0; i < NumberOfParticles; i++) {
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
    int    i, j;
    double xij, yij, zij;
    double distance, distance2;
    double w;

    for (i = 0; i < NumberOfParticles; i++) {
        ParticleNumberDensity[i] = 0.0;
        if (ParticleType[i] == GHOST) continue;
        for (j = 0; j < NumberOfParticles; j++) {
            if ((j == i) || (ParticleType[j] == GHOST)) continue;
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
    int i;
    double n0 = N0_forParticleNumberDensity;
    double beta = THRESHOLD_RATIO_OF_NUMBER_DENSITY;

    for (i = 0; i < NumberOfParticles; i++) {
        //int temp = ParticleType[i];
        if (ParticleType[i] == GHOST || ParticleType[i] == DUMMY_WALL) {
            BoundaryCondition[i] = GHOST_OR_DUMMY;
        }
        else if (ParticleNumberDensity[i] < beta * n0) {
            BoundaryCondition[i] = SURFACE_PARTICLE;
            //ParticleType[i] = SURFACE;
        }
        else {
            //ParticleType[i] = temp;
            BoundaryCondition[i] = INNER_PARTICLE;
        }
    }
}


void setSourceTerm(void) {
    int i;
    double n0 = N0_forParticleNumberDensity;
    double gamma = RELAXATION_COEFFICIENT_FOR_PRESSURE;

    for (i = 0; i < NumberOfParticles; i++) {
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
    double xij, yij, zij;
    double distance, distance2;
    double coefficientIJ;
    double n0 = N0_forLaplacian;
    int    i, j;
    double a;
    int n = NumberOfParticles;

    for (i = 0; i < NumberOfParticles; i++) {
        for (j = 0; j < NumberOfParticles; j++) {
            CoefficientMatrix[i * n + j] = 0.0;
        }
    }

    a = 2.0 * DIM / (n0 * Lambda);
    for (i = 0; i < NumberOfParticles; i++) {
        if (BoundaryCondition[i] != INNER_PARTICLE) continue;
        for (j = 0; j < NumberOfParticles; j++) {
            if ((j == i) || (BoundaryCondition[j] == GHOST_OR_DUMMY)) continue;
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
    exceptionalProcessingForBoundaryCondition();
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
    int    i, j, k;
    double c;
    double sumOfTerms;
    int    n = NumberOfParticles;

    for (i = 0; i < n; i++) {
        Pressure[i] = 0.0;
    }
    for (i = 0; i < n - 1; i++) {
        if (BoundaryCondition[i] != INNER_PARTICLE) continue;
        for (j = i + 1; j < n; j++) {
            if (BoundaryCondition[j] == GHOST_OR_DUMMY) continue;
            c = CoefficientMatrix[j * n + i] / CoefficientMatrix[i * n + i];
            for (k = i + 1; k < n; k++) {
                CoefficientMatrix[j * n + k] -= c * CoefficientMatrix[i * n + k];
            }
            SourceTerm[j] -= c * SourceTerm[i];
        }
    }
    for (i = n - 1; i >= 0; i--) {
        if (BoundaryCondition[i] != INNER_PARTICLE) continue;
        sumOfTerms = 0.0;
        for (j = i + 1; j < n; j++) {
            if (BoundaryCondition[j] == GHOST_OR_DUMMY) continue;
            sumOfTerms += CoefficientMatrix[i * n + j] * Pressure[j];
        }
        Pressure[i] = (SourceTerm[i] - sumOfTerms) / CoefficientMatrix[i * n + i];
    }
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


void writeData_inProfFormat(void) {
    int i;
    FILE* fp;
    char fileName[256];

    sprintf(fileName, "output_%04d.prof", FileNumber);
    fp = fopen(fileName, "w");
    fprintf(fp, "%lf\n", Time);
    fprintf(fp, "%d\n", NumberOfParticles);
    for (i = 0; i < NumberOfParticles; i++) {
        fprintf(fp, "%d %lf %lf %lf %lf %lf %lf %lf %lf\n"
            , ParticleType[i], PositionX[i], PositionY[i], PositionZ[i]
            , VelocityX[i], VelocityY[i], VelocityZ[i], Pressure[i], ParticleNumberDensity[i]);
    }
    fclose(fp);
    FileNumber++;
}


void writeData_inVtuFormat(void) {
    int i;
    double absoluteValueOfVelocity;
    FILE* fp;
    char fileName[1024];

    sprintf(fileName, "particle_%04d.vtu", FileNumber);
    fp = fopen(fileName, "w");
    fprintf(fp, "<?xml version='1.0' encoding='UTF-8'?>\n");
    fprintf(fp, "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
    fprintf(fp, "<UnstructuredGrid>\n");
    fprintf(fp, "<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n", NumberOfParticles, NumberOfParticles);
    fprintf(fp, "<Points>\n");
    fprintf(fp, "<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
    for (i = 0; i < NumberOfParticles; i++) {
        fprintf(fp, "%lf %lf %lf\n", PositionX[i], PositionY[i], PositionZ[i]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Points>\n");
    fprintf(fp, "<PointData>\n");
    fprintf(fp, "<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
    for (i = 0; i < NumberOfParticles; i++) {
        fprintf(fp, "%d\n", ParticleType[i]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
    for (i = 0; i < NumberOfParticles; i++) {
        absoluteValueOfVelocity =
            sqrt(VelocityX[i] * VelocityX[i] + VelocityY[i] * VelocityY[i] + VelocityZ[i] * VelocityZ[i]);
        fprintf(fp, "%f\n", (float)absoluteValueOfVelocity);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='Pressure' format='ascii'>\n");
    for (i = 0; i < NumberOfParticles; i++) {
        fprintf(fp, "%f\n", (float)Pressure[i]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray NumberOfComponents='1' type='Float32' Name='ParticleNumberDensity' format='ascii'>\n");
    for (i = 0; i < NumberOfParticles; i++) {
        fprintf(fp, "%f\n", (float)ParticleNumberDensity[i]);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "<Cells>\n");
    fprintf(fp, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
    for (i = 0; i < NumberOfParticles; i++) {
        fprintf(fp, "%d\n", i);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");
    for (i = 0; i < NumberOfParticles; i++) {
        fprintf(fp, "%d\n", i + 1);
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "<DataArray type='UInt8' Name='types' format='ascii'>\n");
    for (i = 0; i < NumberOfParticles; i++) {
        fprintf(fp, "1\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</Cells>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
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
    //return 0;



	//初始化MPS工具
	//InitMPSTool();
	//ParticlesArrange();
	//SetInitialN0();
}

//void MPSWaterParticleGroup::ParticlesArrange()
//{
//	vec3 initPos = transformation.position;
//	vec3 offset = vec3(l0);									//粒子空间分布布局
//	int width, height, depth;								//粒子的长宽高排布
//	width = height = depth = 3;
//
//	height = 4;
//
//	float highest = 0;										//找到最高点位置
//	if (offset.y > 0)
//		highest = initPos.y + (height - 1) * offset.y;
//	else
//		highest = initPos.y;
//	vector<vec3> positions;
//	int n = 0;
//	n = CubeDistribute(positions, initPos, offset, width, height, depth);
//	for (int i = 0; i < n; i++)
//	{
//		//初始化当前粒子并push进去
//		MPSWaterParticle R;
//		//此处指定粒子属性 位置
//		R.position = positions[i];
//		//如果这个粒子在最上面则标记为表面粒子
//		if (positions[i].y == highest)
//			R.particleType = SURFACE;
//		else
//			R.particleType = FLUID;
//		R.index = i;
//		particles.push_back(R);
//	}
//	particleNumber = n;
//
//	//还需要初始化粒子墙以及dummy wall	(可以看成再创建1+2层空心无盖的cube)   ps:要确保容器包围住所有的粒子
//	//wall
//	int containerWidth = width + 2;
//	int containerHeight = height + 3;
//	int containerDepth = depth + 2;
//	initPos -= offset;
//	positions.clear();
//	n = UncoverHollowCubeDistribute(positions, initPos, offset, containerWidth, containerHeight, containerDepth);
//	for (int i = 0; i < n; i++)
//	{
//		MPSWaterParticle R;
//		R.position = positions[i];
//		R.particleType = WALL;
//		R.index = i + particleNumber;				//加上之前已经有的索引偏移量
//		particles.push_back(R);
//	}
//	particleWallNum = n;
//
//	containerWidth += 2;
//	containerHeight += 2;
//	containerDepth += 2;
//	initPos -= offset;
//	positions.clear();
//	n = UncoverHollowCubeDistribute(positions, initPos, offset, containerWidth, containerHeight, containerDepth);
//	for (int i = 0; i < n; i++)
//	{
//		MPSWaterParticle R;
//		R.position = positions[i];
//		R.particleType = WALL;
//		R.index = i + particleNumber;				//加上之前已经有的索引偏移量
//		particles.push_back(R);
//	}
//	particleWallNum += n;
//
//	//dummy
//	int dummyWidth = containerWidth + 2;
//	int dummyHeight = containerHeight + 2;
//	int dummyDepth = containerDepth + 2;
//	positions.clear();
//	initPos -= offset;
//
//	n = UncoverHollowCubeDistribute(positions, initPos, offset, dummyWidth, dummyHeight, dummyDepth);
//	//第二层
//	dummyHeight += 1;
//	dummyWidth += 2;
//	dummyDepth += 2;
//	initPos -= offset;
//
//	n += UncoverHollowCubeDistribute(positions, initPos, offset, dummyWidth, dummyHeight, dummyDepth);
//
//	//dummyHeight += 1;
//	//dummyWidth += 2;
//	//dummyDepth += 2;
//	//initPos -= offset;
//
//	//n += UncoverHollowCubeDistribute(positions, initPos, offset, dummyWidth, dummyHeight, dummyDepth);
//	for (int i = 0; i < n; i++)
//	{
//		MPSWaterParticle R;
//		R.position = positions[i];
//		R.particleType = DUMMY;
//		R.index = i + particleNumber + particleWallNum;				//加上之前已经有的索引偏移量
//		particles.push_back(R);
//	}
//
//	particleDummyNum = n;
//
//	particleTotalNum = particleNumber + particleDummyNum + particleWallNum;
//}

void MPSWaterParticleGroup::InitMPSTool()
{
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
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
	//MarchingCube::GetInstance()->GetMeshData(posArray, verticesInfo, provideNormal, provideTex);

	shaderData->InitVertexBuffer(verticesInfo, provideNormal, provideTex);
}

// 计算对应的3个n0值，计算方法为假设出一个被完全包围在中间的粒子（粒子影响范围内充满粒子），用这个粒子的密度作为初始密度
void MPSWaterParticleGroup::SetInitialN0()
{
	MPSToolFun* tool = MPSToolFun::GetMPSTool();

	vec3 ri(0, 0, 0), rj(0, 0, 0);
	vec3 r_min(-4);
	vec3 r_max(5);
	double dis;

	if (DIMENSION == 2) {
		r_min.z = 0;
		r_max.z = 1;
	}

	//r_min.z = 0;
	//r_max.z = 1;

	for (int x = r_min.x; x < r_max.x; ++x)
	{
		for (int y = r_min.y; y < r_max.y; ++y)
		{
			for (int z = r_min.z; z < r_max.z; ++z)
			{
				if (((x == 0) && (y == 0)) && (z == 0))continue;
				rj = (float)l0 * vec3(x, y, z);
				dis = distance(rj, ri);

				n0ForNumberDensity += tool->WeightFun(dis, RADIUS_FOR_NUMBER_DENSITY);
				n0ForGradient += tool->WeightFun(dis, RADIUS_FOR_GRADIENT);
				//n0ForLambda += dis * dis * tool->WeightFun(dis, RADIUS_FOR_LAPLACIAN);
				n0ForLambda += tool->WeightFun(dis, RADIUS_FOR_LAPLACIAN);

				lambda0 += dis * dis * tool->WeightFun(dis, RADIUS_FOR_LAPLACIAN);
			}
		}
	}

	lambda0 /= n0ForLambda;

	tool->testLambda0 = lambda0;


}

//void MPSWaterParticleGroup::UpdateAdjoin(float range)
//{
//	int num1 = particleNumber + particleWallNum;
//	int num2 = particles.size();
//	//邻接粒子中普通粒子不需要考虑dummy，wall粒子都要考虑
//	for (int i = 0; i < num1; i++)
//	{
//		particles[i].adjoinParticleIndex.clear();
//		if (i < particleNumber)
//		{
//			//普通粒子的处理
//			for (int j = 0; j < num1; j++)
//			{
//				if (j != i)
//				{
//					vec3 pos1 = particles[j].position;
//					if (distance(pos1, particles[i].position) <= range)
//						particles[i].adjoinParticleIndex.push_back(j);
//				}
//			}
//		}
//		else
//		{
//			//wall粒子的处理
//			for (int j = 0; j < num2; j++)
//			{
//				if (j != i)
//				{
//					vec3 pos1 = particles[j].position;
//					if (distance(pos1, particles[i].position) <= range)
//						particles[i].adjoinParticleIndex.push_back(j);
//				}
//			}
//		}
//	}
//
//}









void MPSWaterParticleGroup::Update(float dt)
{
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


    //iTimeStep++;
    //Time += DT;
    //if ((iTimeStep % OUTPUT_INTERVAL) == 0) {
    //    printf("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticles: %d\n", iTimeStep, Time, NumberOfParticles);
    //    writeData_inVtuFormat();
    //    writeData_inProfFormat();
    //}
    //if (Time >= FINISH_TIME) { break; }
    InitBufferData();
}















//void MPSWaterParticleGroup::Update(float dt)
//{
//	shaderData->UpdateMatrix(transformation);
//	//return;
//	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
//	//计算显示步骤（忽略压力项的第一次速度坐标校正）
//
//	vector<vec3> initUArrayWallAndFluid;			// 时间步n的速度集合
//	vector<vec3> initPArrayWallAndFluid;			// 时间步n的位置集合
//	vector<vec3> middleUArray;						//
//	vector<vec3> middlePArray;						//
//
//	vector<vec3> middlePAllArray;					// 所有粒子的中间值坐标（用来计算Wall的中间密度）
//
//	vector<double> Right;					//存储每一个粒子的右端项
//	vector<bool> ignoreArray;				// 压力泊松方程中需要忽略的项（surface）
//
//	for (int i = 0; i < particleNumber + particleWallNum; i++)
//	{
//		initUArrayWallAndFluid.push_back(particles[i].speed);
//		initPArrayWallAndFluid.push_back(particles[i].position);
//	}
//
//	// 计算每个流体粒子的中间速度
//	for (int i = 0; i < particleNumber + particleWallNum; i++)
//	{
//		vec3 tempU = vec3(particles[i].speed);
//		vec3 tempPos = vec3(particles[i].position);
//		if (i < particleNumber)
//		{
//			tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(initUArrayWallAndFluid, initPArrayWallAndFluid, i, n0ForLambda), particles[i].speed);
//			//tempPos = particles[i].position + (float)0.5 * (tempU + particles[i].speed) * (float)DELTA_TIME;
//			tempPos = particles[i].position + tempU * (float)DELTA_TIME;
//		}
//
//
//		middleUArray.push_back(tempU);
//		middlePArray.push_back(tempPos);
//
//		middlePAllArray.push_back(tempPos);
//	}
//
//	// 将dummy补齐
//	for (int i = particleNumber + particleWallNum; i < particleTotalNum; i++)
//		middlePAllArray.push_back(particles[i].position);
//
//	// Fluid的中间粒子数密度
//	for (int i = 0; i < particleNumber; i++)
//	{
//		float tempN = mpsTool->DensityN(middlePArray, i);
//		particles[i].middleN = tempN;
//		//cout << n0ForNumberDensity << " " << tempN << endl;
//		// 表面检测
//		if (tempN < n0ForNumberDensity * THRESHOLD_RATIO_OF_NUMBER_DENSITY)
//		{
//			particles[i].particleType = SURFACE;
//			ignoreArray.push_back(true);
//		}
//		else
//		{
//			particles[i].particleType == FLUID;
//			ignoreArray.push_back(false);
//		}
//
//
//
//
//	}
//	// Wall的中间粒子数密度
//	for (int i = particleNumber; i < particleNumber + particleWallNum; i++)
//	{
//		float tempN = mpsTool->DensityN(middlePAllArray, i);
//		particles[i].middleN = tempN;
//		ignoreArray.push_back(false);
//		if (tempN < n0ForNumberDensity * THRESHOLD_RATIO_OF_NUMBER_DENSITY)
//		{
//			particles[i].particleType = SURFACE;
//			ignoreArray.push_back(true);
//		}
//		else
//		{
//			particles[i].particleType == WALL;
//			ignoreArray.push_back(false);
//		}
//	}
//
//
//	//for (int i = 0; i < particleNumber + particleWallNum; i++)
//	//{
//	//	particles[i].speed = middleUArray[i];
//	//	particles[i].position = middlePArray[i];
//	//}
//
//	//InitBufferData();
//	//return;
//	vector<double> SourceTerm;
//	SourceTerm.resize(5000, 0);
//
//	int i;
//	double n0 = n0ForNumberDensity;
//	double NumberOfParticles = particles.size();
//	double gamma = 0.2;
//
//	for (i = 0; i < NumberOfParticles; i++) {
//		SourceTerm[i] = 0.0;
//		if (particles[i].particleType == DUMMY) continue;
//		if (particles[i].particleType == WALL || particles[i].particleType == FLUID) {
//			SourceTerm[i] = gamma * (1.0 / (DELTA_TIME * DELTA_TIME)) * ((particles[i].middleN - n0) / n0);
//		}
//		else if (particles[i].particleType == SURFACE) {
//			SourceTerm[i] = 0.0;
//		}
//	}
//
//
//
//	vector<double> CoefficientMatrix;
//
//
//	CoefficientMatrix.resize(25000000, 0);
//
//
//	double xij, yij, zij;
//	double distance, distance2;
//	double coefficientIJ;
//
//
//	int j;
//	double a;
//	int n = particles.size();
//
//	for (i = 0; i < NumberOfParticles; i++) {
//		for (j = 0; j < NumberOfParticles; j++) {
//			CoefficientMatrix[i * n + j] = 0.0;
//		}
//	}
//
//	a = 2.0 * DIMENSION / (n0 * lambda0);
//	for (i = 0; i < NumberOfParticles; i++) {
//		if (particles[i].particleType == WALL || particles[i].particleType == FLUID)
//		{
//			for (j = 0; j < NumberOfParticles; j++) {
//				if ((j == i) || (particles[j].particleType == DUMMY)) continue;
//				xij = particles[j].position.x - particles[i].position.x;
//				yij = particles[j].position.y - particles[i].position.y;
//				zij = particles[j].position.z - particles[i].position.z;
//				distance2 = (xij * xij) + (yij * yij) + (zij * zij);
//				distance = sqrt(distance2);
//				if (distance >= RADIUS_FOR_LAPLACIAN)continue;
//				coefficientIJ = a * mpsTool->WeightFun(distance, RADIUS_FOR_LAPLACIAN) / FLUID_DENSITY;
//				CoefficientMatrix[i * n + j] = (-1.0) * coefficientIJ;
//				CoefficientMatrix[i * n + i] += coefficientIJ;
//			}
//			CoefficientMatrix[i * n + i] += (0.45E-9) / (DELTA_TIME * DELTA_TIME);
//
//		}
//		else
//		{
//			continue;
//		}
//;
//	}
//
//	int k;
//	double c;
//	double sumOfTerms;
//	//int    n = NumberOfParticles;
//
//	for (i = 0; i < n; i++) {
//		particles[i].pressure = 0.0;
//	}
//	for (i = 0; i < n - 1; i++) {
//
//		if (particles[i].particleType == WALL || particles[i].particleType == FLUID)
//		{
//			for (j = i + 1; j < n; j++) {
//				if (particles[j].particleType == DUMMY) continue;
//				c = CoefficientMatrix[j * n + i] / CoefficientMatrix[i * n + i];
//				for (k = i + 1; k < n; k++) {
//					CoefficientMatrix[j * n + k] -= c * CoefficientMatrix[i * n + k];
//				}
//				SourceTerm[j] -= c * SourceTerm[i];
//			}
//		}
//		//if (BoundaryCondition[i] != INNER_PARTICLE) continue;
//	}
//	for (i = n - 1; i >= 0; i--) {
//		if (particles[i].particleType == WALL || particles[i].particleType == FLUID)
//		{
//			sumOfTerms = 0.0;
//			for (j = i + 1; j < n; j++) {
//				if (particles[j].particleType == DUMMY) continue;
//				sumOfTerms += CoefficientMatrix[i * n + j] * particles[j].pressure;
//			}
//			particles[i].pressure = (SourceTerm[i] - sumOfTerms) / CoefficientMatrix[i * n + i];
//		}
//
//	}
//
//	vector<double> resP;
//	for (int i = 0; i < particles.size(); i++)
//	{
//		resP.push_back(particles[i].pressure);
//	}
//
//	for (int i = 0; i < particleNumber + particleWallNum; i++)
//	{
//		vec3 secondTempU = middleUArray[i];
//		vec3 secondTempP = middlePArray[i];
//		if (i < particleNumber)
//		{
//			vec3 resGradient = mpsTool->OldGradient(middlePArray, resP, i, n0ForGradient);
//
//			secondTempU = (float)(-DELTA_TIME / FLUID_DENSITY) * resGradient + middleUArray[i];
//			//secondTempP = middlePArray[i] + (float)0.5 * (middleUArray[i] + secondTempU) * (float)DELTA_TIME;
//			//secondTempP = middlePArray[i] + secondTempU * (float)DELTA_TIME;
//			secondTempP = middlePArray[i] + (secondTempU - middleUArray[i]) * (float)DELTA_TIME;
//		}
//
//		particles[i].speed = secondTempU;
//		particles[i].position = secondTempP;
//
//	}
//
//
//
//
//
//	// Wall和Fluid粒子的右端项（其中要排除对surface的计算）
//	for (int i = 0; i < particleNumber + particleWallNum; i++)
//	{
//		float resRight = 0;
//		//if (!ignoreArray[i])
//		//{
//
//		//}
//		resRight = mpsTool->OldImplicitLaplacianRight(FLUID_DENSITY, n0ForNumberDensity, particles[i].middleN);
//		Right.push_back(resRight);
//	}
//
//
//
//	// 解一个泊松方程
//	//vector<double> resP = mpsTool->ImplicitCalculateP(middlePArray, n0ForLambda, ignoreArray, Right);
//
//	for (int i = 0; i < resP.size(); i++)
//	{
//		if (ignoreArray[i])
//			resP[i] = 0;
//	}
//
//
//	for (int i = 0; i < particleNumber + particleWallNum; i++)
//	{
//		vec3 secondTempU = middleUArray[i];
//		vec3 secondTempP = middlePArray[i];
//		if (i < particleNumber)
//		{
//			vec3 resGradient = mpsTool->OldGradient(middlePArray, resP, i, n0ForGradient);
//
//			secondTempU = (float)(-DELTA_TIME / FLUID_DENSITY) * resGradient + middleUArray[i];
//			//secondTempP = middlePArray[i] + (float)0.5 * (middleUArray[i] + secondTempU) * (float)DELTA_TIME;
//			//secondTempP = middlePArray[i] + secondTempU * (float)DELTA_TIME;
//			secondTempP = middlePArray[i] + (secondTempU - middleUArray[i]) * (float)DELTA_TIME;
//		}
//
//		particles[i].speed = secondTempU;
//		particles[i].position = secondTempP;
//
//	}
//
//	printf("sds");
//
//	//将计算好的位置点进行建模
//	//Modeling();
//	InitBufferData();
//}

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
                color = vec3(0, 255, 0);
            break;
        case SURFACE_PARTICLE:
            color = vec3(255, 0, 0);
            break;
        case GHOST_OR_DUMMY:
            //color = vec3(0);
            color = vec3(0, 255, 255);
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
