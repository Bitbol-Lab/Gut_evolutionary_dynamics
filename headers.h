#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define pi 3.14159265359
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// Structure to hold all simulation parameters

typedef struct 
{
    // Simulation parameters
    double dt;           // Time step
    double dx;           // Spatial step
    long int t_init;     // Number of time steps for initialization 
    int t_inter;     // Number of time steps for initialization between successive stochastic simulation
    int n_dx;           // Maximum spatial position index
    double D_0;          // Diffusion coefficient
    double v_0;          // Drift velocity
    double r_wild;       // Growth rate wild-type
    double r_mutant;     // Growth rate mutant
    double F_in;         // Food influx
    double k_rep;        // Monod constant
    double alpha;        // Food consumption coefficient
    double conv_threshold; // Convergence threshold for deterministic initialization
    int x_init; 
    int n_simu; 
    int simu_id;
    
    // Derived parameters (computed from above)
    double p_diff;
    double p_drift;
    double p_left;
    double p_right;
    double p_right_0;
    double p_left_max;
    double p_right_max;
    double *p_array;
    double *p_array_max;

} SimParams;


void replace_char(char *str, char find, char replace);
double proba_repr(double r_rate, double food, const SimParams *params);
double compute_n_active(int *wild_array, double *food_array, const SimParams *params);
void generate_p_mutant(double *p_repr_mutant, double *food_array, const SimParams *params);

void func_food_transport(double *food_array, double *food_array_next, const SimParams *params);
void update_transport_wild_only(int *n_wild, int *wild_array, int *wild_array_next, const SimParams *params, gsl_rng *r);
void update_replication_wild_only(int *n_wild, int *wild_array, double *food_array, const SimParams *params, gsl_rng *r);
void update_replication_first_mutant(int *n_wild, int *wild_array, int *n_mutant, int *mutant_array, double *food_array, const SimParams *params, gsl_rng *r);
void update_transport_wild_mutant(int *n_wild, int *wild_array, int *wild_array_next, int *n_mutant, int *mutant_array, int *mutant_array_next, const SimParams *params, gsl_rng *r);
void update_replication_wild_mutant(int *n_wild, int *wild_array, int *n_mutant, int *mutant_array, double *food_array, double *food_array_next, const SimParams *params, gsl_rng *r);

double generate_c_wild_food(double *c_wild_array, double *food_array, const SimParams *params);
void generate_array_bacteria(int n_wild, int *wild_array, double *c_wild_array, const SimParams *params, gsl_rng *r);
void main_func_init(const char *nom_fichier1,const char *nom_fichier2, const SimParams *params, gsl_rng *r);
void extract_data(const char *nom_fichier, int *n_wild_init, int *wild_array_init, double *food_array_init, const SimParams *params);
void main_func_phi_moy(const char *name_file, const SimParams *params, gsl_rng *r);
void main_func_phi_x(const char *nom_fichier, const SimParams *params, gsl_rng *r);
