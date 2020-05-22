#pragma once
//constexpr unsigned int DIMENSION = 2;
//constexpr double PARTICLE_DISTANCE = 0.025;
//constexpr double DELTA_TIME = 0.001;

constexpr unsigned int DIMENSION										= 3;
constexpr double PARTICLE_DISTANCE										= 0.05; 
constexpr double DELTA_TIME												= 0.005;

constexpr double FINISH_TIME											= 1.0;
constexpr double KINEMATIC_VISCOSITY									= 1.0E-6;					// 粘度系数
constexpr double FLUID_DENSITY											= 1000.0;
constexpr double GRAVITY_X												= 0.0;
constexpr double GRAVITY_Y												= -9.8;
constexpr double GRAVITY_Z												= 0.0;
constexpr double RADIUS_FOR_NUMBER_DENSITY								= 2.1 * PARTICLE_DISTANCE;
constexpr double RADIUS_FOR_GRADIENT									= 2.1 * PARTICLE_DISTANCE;
constexpr double RADIUS_FOR_LAPLACIAN									= 3.1 * PARTICLE_DISTANCE;
constexpr double COLLISION_DISTANCE										= 0.5 * PARTICLE_DISTANCE;
constexpr double THRESHOLD_RATIO_OF_NUMBER_DENSITY						= 0.97;						// 表面检测的beta值
constexpr double RELAXATION_COEFFICIENT_FOR_PRESSURE					= 0.2;						// gamma值


//constexpr double COEFFICIENT_OF_RESTITUTION = 0.2; // 斀敪學悢
//constexpr double COMPRESSIBILITY = 0.45E-9;
//constexpr double EPS = 0.01 * PARTICLE_DISTANCE;