#pragma once
//constexpr unsigned int DIMENSION = 2;
//constexpr double PARTICLE_DISTANCE = 0.025;
//constexpr double DELTA_TIME = 0.001;

//constexpr unsigned int DIMENSION										= 3;
//constexpr double PARTICLE_DISTANCE										= 0.05; 
//constexpr double DELTA_TIME												= 0.005;
//
//constexpr double FINISH_TIME											= 1.0;
//constexpr double KINEMATIC_VISCOSITY									= 1.0E-6;					// 粘度系数
//constexpr double FLUID_DENSITY											= 1000.0;
//constexpr double GRAVITY_X												= 0.0;
//constexpr double GRAVITY_Y												= -9.8;
//constexpr double GRAVITY_Z												= 0.0;
//constexpr double RADIUS_FOR_NUMBER_DENSITY								= 2.1 * PARTICLE_DISTANCE;
//constexpr double RADIUS_FOR_GRADIENT									= 2.1 * PARTICLE_DISTANCE;
//constexpr double RADIUS_FOR_LAPLACIAN									= 3.1 * PARTICLE_DISTANCE;
//constexpr double COLLISION_DISTANCE										= 0.5 * PARTICLE_DISTANCE;
//constexpr double THRESHOLD_RATIO_OF_NUMBER_DENSITY						= 0.97;						// 表面检测的beta值
//constexpr double RELAXATION_COEFFICIENT_FOR_PRESSURE					= 0.2;						// gamma值


//constexpr double COEFFICIENT_OF_RESTITUTION = 0.2; // 斀敪學悢
//constexpr double COMPRESSIBILITY = 0.45E-9;
//constexpr double EPS = 0.01 * PARTICLE_DISTANCE;


constexpr unsigned int DIM								= 2;
constexpr float PARTICLE_DISTANCE						= 0.01;
constexpr float DT										= 0.002;
constexpr float OUTPUT_INTERVAL							= 20;

//constexpr float DIM									= 3;
//constexpr float PARTICLE_DISTANCE						= 0.075;
//constexpr float DT									= 0.003;
//constexpr float OUTPUT_INTERVAL						= 2;

constexpr unsigned int ARRAY_SIZE						= 5000;
//#define FINISH_TIME          2.0
constexpr float KINEMATIC_VISCOSITY						= 1.0E-6;
constexpr float FLUID_DENSITY							= 1000.0;
constexpr float GRAVITY_X								= 0.0;
constexpr float GRAVITY_Y								= -9.8;
constexpr float GRAVITY_Z								= 0.0;
constexpr float RADIUS_FOR_NUMBER_DENSITY				= 2.1 * PARTICLE_DISTANCE;
constexpr float RADIUS_FOR_GRADIENT						= 2.1 * PARTICLE_DISTANCE;
constexpr float RADIUS_FOR_LAPLACIAN					= 3.1 * PARTICLE_DISTANCE;
constexpr float COLLISION_DISTANCE						= 0.5 * PARTICLE_DISTANCE;
constexpr float THRESHOLD_RATIO_OF_NUMBER_DENSITY		= 0.97;
constexpr float COEFFICIENT_OF_RESTITUTION				= 0.2;
constexpr float COMPRESSIBILITY							= 0.45E-9;
constexpr float EPS										= 0.01 * PARTICLE_DISTANCE;
constexpr float RELAXATION_COEFFICIENT_FOR_PRESSURE = 0.2;
constexpr int ON										= 1;
constexpr int OFF										= 0;

constexpr int GHOST										= -1;
constexpr int FLUID										= 0;
constexpr int SURFACE									= 1;
constexpr int WALL										= 2;
constexpr int DUMMY_WALL								= 3;
constexpr int GHOST_OR_DUMMY							= -1;
constexpr int SURFACE_PARTICLE							= 1;
constexpr int INNER_PARTICLE							= 0;
constexpr int DIRICHLET_BOUNDARY_IS_NOT_CONNECTED		= 0;
constexpr int DIRICHLET_BOUNDARY_IS_CONNECTED			= 1;
constexpr int DIRICHLET_BOUNDARY_IS_CHECKED				= 2;
