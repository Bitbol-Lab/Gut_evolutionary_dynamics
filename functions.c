#include "headers.h"

void replace_char(char *str, char find, char replace) 
{
    /* Replace every occurrence of a specified character in a string. */ 

    for (char *p = str; *p; p++) 
    {
        if (*p == find) 
        {
            *p = replace;
        }
    }
}

double proba_repr(double r_rate, double food, const SimParams *params) 
{
    /* Compute the Michaelis-Menten replication probability over one time step.*/
     if (food <= 0) 
    {
        return 0.0; 
    }
    double p= r_rate * params->dt * food / (params->k_rep + food);
    if (p> 1.)
    {
        return 1.;
    }
    return p;
}

double compute_n_active(int *wild_array, double *food_array, const SimParams *params)
{
    /* Compute the active population size, as in Ref.  
    Labavić, D., Loverdo, C., & Bitbol, A. F. (2022). Hydrodynamic flow and concentration gradients in the gut enhance neutral bacterial diversity. Proceedings of the National Academy of Sciences, 119(1), e2108671119. */ 
    double n_active = 0;
    double rho_0 = (food_array[0]) / (food_array[0]+ params->k_rep)*0.5;

    for (int k= 0; k<= params->n_dx; k++) 
    {
        if ((food_array[k]) / (food_array[k]+ params->k_rep) >= rho_0)
        {
            n_active += wild_array[k] ;
        }
    }
    return n_active;
}

void func_food_transport(double *food_array, double *food_array_next, const SimParams *params) 
{
    /*Update food concentrations by advection and diffusion. */
    for (int k = 1; k < params->n_dx; k++) 
    {
        food_array_next[k]= food_array[k] + params->p_diff * ( food_array[k+1] - 2 * food_array[k]  + food_array[k-1] ) - params->p_drift /2 * ( food_array[k+1] - food_array[k-1]);
    }
    
    // Apply inlet and outlet boundary conditions.
    food_array_next[0] = food_array[0] + 2 * params->p_diff * ( food_array[1] - food_array[0] ) - (2* params->p_drift + pow(params->p_drift, 2) / params->p_diff)* (food_array[0] - params->F_in);
    food_array_next[params->n_dx] = food_array[params->n_dx] + 2 * params->p_diff * ( food_array[params->n_dx-1] - food_array[params->n_dx]);

    // Store the updated food profile.
    memcpy(food_array, food_array_next, (params->n_dx + 1) * sizeof(double));
}

void update_transport_wild_only(int *n_wild, int *wild_array, int *wild_array_next, const SimParams *params, gsl_rng *r)
{
    /*Stochastically transport wild-type bacteria, accounting for downstream washout. */
    unsigned int count_n[3];
    int n, n_right;

    memset(wild_array_next, 0, (params->n_dx + 1)*sizeof(int));

    // Sample transport at interior sites.
    for (int k= 1; k< params->n_dx; k++)
    {   
        n = wild_array[k];
        if (n>0) 
        {
            gsl_ran_multinomial(r, 3, n, params->p_array, count_n);
            wild_array_next[k] += count_n[0];
            wild_array_next[k-1] += count_n[1];
            wild_array_next[k+1] += count_n[2];  
        }
    }

    // Reflect bacteria at x = 0.
    n = wild_array[0];
    if (n>0)
    {
        n_right = gsl_ran_binomial(r, params->p_right_0, n);
        wild_array_next[0] += n- n_right;
        wild_array_next[1] += n_right;
    }

    // Allow bacteria to leave the system at x = L.  
    n = wild_array[params->n_dx];
    if (n>0)
    {
        gsl_ran_multinomial(r, 3, n, params->p_array_max, count_n);
        wild_array_next[params->n_dx] += count_n[0];
        wild_array_next[params->n_dx-1] += count_n[1];
        (*n_wild)= (*n_wild)- count_n[2]; 
    }

    memcpy(wild_array, wild_array_next, (params->n_dx + 1) * sizeof(int));
}


void update_replication_wild_only(int *n_wild, int *wild_array, double *food_array, const SimParams *params, gsl_rng *r)
{
    /* Sample wild-type divisions and subtract the corresponding food consumption.*/
    int n, n_repr;
    double p_repr;


    for (int k= 0; k<= params->n_dx; k++)
    {          
        n = wild_array[k];
        if (n>0) 
        {
            p_repr = proba_repr(params->r_wild, food_array[k], params);
            n_repr = gsl_ran_binomial(r, p_repr, n);
            (*n_wild)= (*n_wild) +n_repr; 
            wild_array[k] = wild_array[k] + n_repr;
            food_array[k] = food_array[k] - n_repr *1.0 / params->alpha;
        }
    }
}


void update_replication_first_mutant(int *n_wild, int *wild_array, int *n_mutant, int *mutant_array, double *food_array, const SimParams *params, gsl_rng *r)
{
    /* Sample wild-type divisions and introduce one mutant among the new offspring.
    If several divisions occur, select one uniformly among these events.  */
    int n, n_repr;
    double p_repr;
    int n_repr_tot = 0;
 
    double* distr_replication = calloc((params->n_dx + 1), sizeof(double));

    // Sample replication events at each site.
    for (int k= 0; k<= params->n_dx; k++)
    {          
        n = wild_array[k];
        if (n>0) 
        {
            p_repr = proba_repr(params->r_wild, food_array[k], params);
            n_repr = gsl_ran_binomial(r, p_repr, n);
            (*n_wild)= (*n_wild) +n_repr;
            wild_array[k] = wild_array[k] + n_repr;
            food_array[k] = food_array[k] - n_repr *1.0 / params->alpha;
            distr_replication[k] += n_repr;
            n_repr_tot += n_repr; 
        }
    }

    // Select one newborn uniformly among all replication events.
    if (n_repr_tot> 0)
    {
        gsl_ran_discrete_t *distr= gsl_ran_discrete_preproc(params->n_dx+1, distr_replication);
        int pos_mutant = (int)gsl_ran_discrete(r, distr);
        
        (*n_mutant)= 1;
        mutant_array[pos_mutant]++;

        (*n_wild)--;
        wild_array[pos_mutant]--;

        gsl_ran_discrete_free(distr);
    }
    free(distr_replication);
}


void update_transport_wild_mutant(int *n_wild, int *wild_array, int *wild_array_next, int *n_mutant, int *mutant_array, int *mutant_array_next, const SimParams *params, gsl_rng *r)
{
    /* Transport wild-type and mutant bacteria independently, including downstream washout. */
    unsigned int count_n[3];
    int n, n_right;


    memset(wild_array_next, 0, (params->n_dx + 1)*sizeof(int));
    memset(mutant_array_next, 0, (params->n_dx + 1)*sizeof(int));

    // Sample transport at interior sites.
    for (int k= 1; k< params->n_dx; k++)
    {       
        n = wild_array[k];
        if (n>0) 
        {
            gsl_ran_multinomial(r, 3, n, params->p_array, count_n);
            wild_array_next[k] += count_n[0];
            wild_array_next[k-1] += count_n[1];
            wild_array_next[k+1] += count_n[2];  
        }
        n = mutant_array[k];
        if (n>0) 
        {
            gsl_ran_multinomial(r, 3, n, params->p_array, count_n);
            mutant_array_next[k] += count_n[0];
            mutant_array_next[k-1] += count_n[1];
            mutant_array_next[k+1] += count_n[2];  
        }
    }

    // Reflect both populations at x = 0.
    n = wild_array[0];
    if (n>0)
    {
        n_right = gsl_ran_binomial(r, params->p_right_0, n);
        wild_array_next[0] += n- n_right;
        wild_array_next[1] += n_right;
    }
    n = mutant_array[0];
    if (n>0)
    {
        n_right = gsl_ran_binomial(r, params->p_right_0, n);
        mutant_array_next[0] += n- n_right;
        mutant_array_next[1] += n_right;
    }

    // Account for washout of both populations at the outlet x = L.
    n = wild_array[params->n_dx];
    if (n>0)
    {
        gsl_ran_multinomial(r, 3, n, params->p_array_max, count_n);
        wild_array_next[params->n_dx] += count_n[0];
        wild_array_next[params->n_dx-1] += count_n[1];
        (*n_wild)= (*n_wild)- count_n[2]; 
    }
    n = mutant_array[params->n_dx];
    if (n>0)
    {
        gsl_ran_multinomial(r, 3, n, params->p_array_max, count_n);
        mutant_array_next[params->n_dx] += count_n[0];
        mutant_array_next[params->n_dx-1] += count_n[1];
        (*n_mutant)= (*n_mutant)- count_n[2]; 
    }
    memcpy(wild_array, wild_array_next, (params->n_dx + 1) * sizeof(int));
    memcpy(mutant_array, mutant_array_next, (params->n_dx + 1) * sizeof(int));
}

void update_replication_wild_mutant(int *n_wild, int *wild_array, int *n_mutant, int *mutant_array, double *food_array, double *food_array_next, const SimParams *params, gsl_rng *r)
{
    /* Sample wild-type and mutant divisions using their respective growth rates and shared food. */
    int n, n_repr;
    double p_repr;

    memcpy(food_array_next, food_array, (params->n_dx + 1)* sizeof(double));


    for (int k= 0; k<= params->n_dx; k++)
    {      
        // Sample wild-type divisions.    
        n = wild_array[k];
        if (n>0) 
        {
            p_repr = proba_repr(params->r_wild, food_array[k], params);
            n_repr = gsl_ran_binomial(r, p_repr, n);
            (*n_wild)= (*n_wild) +n_repr; 
            wild_array[k] = wild_array[k] + n_repr;
            food_array_next[k] = food_array_next[k] - n_repr *1.0 / params->alpha;
        }

        // Sample mutant divisions.    
        n =mutant_array[k];
        if (n>0) 
        {
            p_repr = proba_repr(params->r_mutant, food_array[k], params);
            n_repr = gsl_ran_binomial(r, p_repr, n);
            (*n_mutant)= (*n_mutant) +n_repr; 
            mutant_array[k] = mutant_array[k] + n_repr;
            food_array_next[k] = food_array_next[k] - n_repr *1.0 / params->alpha;
        }
    }    
    memcpy(food_array, food_array_next, (params->n_dx + 1) * sizeof(double));

}


double generate_c_wild_food(double *c_wild_array, double *food_array, const SimParams *params)
{
    /*Relax the deterministic bacterial and food profiles toward a steady state using explicit Euler steps.  */
    double* c_wild_array_next= malloc((params->n_dx+1)* sizeof(double));
    double* food_array_next= malloc((params->n_dx+1)* sizeof(double));
    
    double diffusion_rate = params->D_0 *params->dt / pow(params->dx, 2);
    double advection_rate = params->v_0 *params->dt / params->dx;

    for (int i = 0; i<= params->n_dx; i++)
    {
        c_wild_array[i] = i* 1./ params->n_dx ;
        food_array[i] = params->F_in* (1- i *1./ params->n_dx) ;
    }

    double total_abs_change =  1 +params->conv_threshold;

    while (total_abs_change > params->conv_threshold)
    {   
        for (int k = 1 ; k < params->n_dx; k++)
        {
            c_wild_array_next[k] = c_wild_array[k] + diffusion_rate* (c_wild_array[k+1] - 2 * c_wild_array[k]  + c_wild_array[k-1]) - advection_rate /2 * (c_wild_array[k+1] - c_wild_array[k-1]) + params->r_wild* params->dt* food_array[k]* c_wild_array[k] /(params->k_rep +food_array[k]);
            food_array_next[k] = food_array[k] + diffusion_rate* (food_array[k+1] - 2 * food_array[k]  + food_array[k-1]) - advection_rate /2 * (food_array[k+1] - food_array[k-1]) - params->r_wild* params->dt* 1.0/ params->alpha * food_array[k]* c_wild_array[k] /(params->k_rep + food_array[k]);
        }
        c_wild_array_next[0] = c_wild_array[0] + 2 * diffusion_rate* (c_wild_array[1] - c_wild_array[0]) - (2 * advection_rate + pow( advection_rate, 2) / diffusion_rate)* c_wild_array[0] + params->r_wild* params->dt* food_array[0]* c_wild_array[0]/(params->k_rep + food_array[0]);
        food_array_next[0] = food_array[0] + 2 * diffusion_rate* (food_array[1] - food_array[0]) - (2* advection_rate + pow( advection_rate, 2) / diffusion_rate)* (food_array[0] - params->F_in) - params->r_wild* params->dt* 1.0/ params->alpha* food_array[0]* c_wild_array[0]/(params->k_rep + food_array[0]);

        c_wild_array_next[params->n_dx] = c_wild_array[params->n_dx] + 2 * diffusion_rate* (c_wild_array[params->n_dx-1] - c_wild_array[params->n_dx]) + params->r_wild* params->dt* food_array[params->n_dx]* c_wild_array[params->n_dx] / (params->k_rep +food_array[params->n_dx]);
        food_array_next[params->n_dx] = food_array[params->n_dx] + 2 * diffusion_rate* (food_array[params->n_dx-1] - food_array[params->n_dx]) - params->r_wild* params->dt* 1.0/params->alpha * food_array[params->n_dx]* c_wild_array[params->n_dx] / (params->k_rep +food_array[params->n_dx]);
        
        total_abs_change = 0;
        for (int i = 0; i <= params->n_dx; i++) 
        {
            total_abs_change += fabs(c_wild_array_next[i] - c_wild_array[i]);
            total_abs_change += fabs(food_array_next[i] - food_array[i]);
        }   
        
        for (int i = 0; i<= params->n_dx; i++)
        {
            c_wild_array[i]= c_wild_array_next[i];
            food_array[i]= food_array_next[i];
        }
    }

    double n_wild = 0;
    for (int k = 0; k<= params->n_dx; k++)
    {
        n_wild = n_wild + c_wild_array[k];
    }

    free(c_wild_array_next);
    free(food_array_next);

    return n_wild;
}

void generate_array_bacteria(int n_wild, int *wild_array, double *c_wild_array, const SimParams *params, gsl_rng *r)
{
    /* Sample discrete bacterial positions from the deterministic concentration profile. */
    int i;
    gsl_ran_discrete_t *distr= gsl_ran_discrete_preproc(params->n_dx+1, c_wild_array);

    for (int k = 0; k < n_wild; k++)
    {
        i = (int)gsl_ran_discrete(r, distr);
        wild_array[i]++;
    }
    gsl_ran_discrete_free(distr);
}



void main_func_init(const char *name_file1, const char *name_file2, const SimParams *params, gsl_rng *r)
{
    /* Compute deterministic steady-state profiles, then initialize
       and evolve a stochastic wild-type population. */

    long now1, now2, duration;
    int hours, minutes, seconds;
    now1 = time(NULL);
    FILE* file;

    // Allocate arrays for bacterial populations and food concentrations.
    double* c_wild_array = calloc((params->n_dx+1), sizeof(double));
    int* wild_array = calloc((params->n_dx+1), sizeof(int));
    int* wild_array_next = calloc((params->n_dx+1), sizeof(int));
    double* food_array = calloc((params->n_dx+1), sizeof(double));
    double* food_array_next = calloc((params->n_dx+1), sizeof(double));

    // Compute deterministic steady-state profiles.
    int n_wild = generate_c_wild_food(c_wild_array, food_array, params);

    now2 = time(NULL);
    duration = now2 - now1;

    hours = duration / 3600;
    minutes = (duration - hours * 3600) / 60;
    seconds = (duration - hours * 3600) % 60;

    // Save deterministic profiles and simulation parameters.
    file = fopen(name_file1, "w");

    fprintf(file, "Time: %d hours, %d minutes, %d seconds. \n\n", hours, minutes, seconds);

    fprintf(file, "dt\t%e\n", params->dt);
    fprintf(file, "dx\t%e\n", params->dx);
    fprintf(file, "t_init\t%ld\n", params->t_init);
    fprintf(file, "n_dx\t%d\n", params->n_dx);
    fprintf(file, "D_0\t%e\n", params->D_0);
    fprintf(file, "v_0\t%e\n", params->v_0);
    fprintf(file, "r_wild\t%e\n", params->r_wild);
    fprintf(file, "F_in\t%e\n", params->F_in);
    fprintf(file, "k_rep\t%e\n", params->k_rep);
    fprintf(file, "alpha\t%e\n", params->alpha);
    fprintf(file, "conv_threshold\t%e\n", params->conv_threshold);
    fprintf(file, "n_wild\t%d\n\n", n_wild);
    fprintf(file, "Cwild_array\tfood_array\n\n");

    for (int k = 0; k <= params->n_dx; k++)
    {
        fprintf(file, "%e\t%e\n", c_wild_array[k], food_array[k]);
    }
    fclose(file);

    // Sample discrete bacterial positions from the steady-state profile.
    generate_array_bacteria(n_wild, wild_array, c_wild_array, params, r);

    // Evolve the stochastic population to relax toward stationarity.
    long int t = 0;
    while (t < params->t_init)
    {
        update_transport_wild_only(&n_wild, wild_array, wild_array_next, params, r);
        update_replication_wild_only(&n_wild, wild_array, food_array, params, r);
        func_food_transport(food_array, food_array_next, params);

        t++;
        if (n_wild == 0)
        {
            fprintf(stderr, "Warning: Wash-out!\n");
        }
    }

    now2 = time(NULL);
    duration = now2 - now1;

    hours = duration / 3600;
    minutes = (duration - hours * 3600) / 60;
    seconds = (duration - hours * 3600) % 60;

    // Save the stochastic initial configuration for subsequent simulations.
    file = fopen(name_file2, "w");

    fprintf(file, "Time: %d hours, %d minutes, %d seconds. \n\n", hours, minutes, seconds);

    fprintf(file, "dt\t%e\n", params->dt);
    fprintf(file, "dx\t%e\n", params->dx);
    fprintf(file, "t_init\t%ld\n", params->t_init);
    fprintf(file, "n_dx\t%d\n", params->n_dx);
    fprintf(file, "D_0\t%e\n", params->D_0);
    fprintf(file, "v_0\t%e\n", params->v_0);
    fprintf(file, "r_wild\t%e\n", params->r_wild);
    fprintf(file, "F_in\t%e\n", params->F_in);
    fprintf(file, "k_rep\t%e\n", params->k_rep);
    fprintf(file, "alpha\t%e\n", params->alpha);
    fprintf(file, "conv_threshold\t%e\n", params->conv_threshold);
    fprintf(file, "n_wild\t%d\n\n", n_wild);
    fprintf(file, "n_wild_array\tfood_array\n\n");

    for (int k = 0; k <= params->n_dx; k++)
    {
        fprintf(file, "%d\t%e\n", wild_array[k], food_array[k]);
    }

    fclose(file);

    printf("\n Time: %d hours, %d minutes, %d seconds. \n", hours, minutes, seconds);

    // Free allocated memory.
    free(c_wild_array);
    free(wild_array);
    free(wild_array_next);
    free(food_array);
    free(food_array_next);
}




void extract_data(const char *name_file, int *n_wild_init, int *wild_array_init, double *food_array_init, const SimParams *params)
{
    /* Load wild-type counts and food concentrations from a file after checking parameter consistency.*/
    FILE *file = fopen(name_file, "r");
    if (file == NULL)
    {
        printf("Error opening file.\n");
        exit(EXIT_FAILURE);
    }

    char buffer[256];
    SimParams params_test;
    int n_wild_value;

    // Skip the runtime header.
    if (!fgets(buffer, sizeof(buffer), file) || !fgets(buffer, sizeof(buffer), file))
    {
        fprintf(stderr, "Error reading header lines\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    char name[64];

    // Parse and validate the stored simulation parameters.
    #define READ_PARAM(fmt, varptr, label) \
        do { \
            if (!fgets(buffer, sizeof(buffer), file)) { \
                fprintf(stderr, "Error reading parameter '%s'\n", label); \
                fclose(file); exit(EXIT_FAILURE); \
            } \
            if (sscanf(buffer, "%63s\t" fmt, name, varptr) != 2 || strcmp(name, label) != 0) { \
                fprintf(stderr, "Error parsing parameter '%s' (line: %s)\n", label, buffer); \
                fclose(file); exit(EXIT_FAILURE); \
            } \
        } while (0)

    READ_PARAM("%lf", &params_test.dt,             "dt");
    READ_PARAM("%lf", &params_test.dx,             "dx");
    READ_PARAM("%ld", &params_test.t_init,         "t_init");
    READ_PARAM("%d",  &params_test.n_dx,           "n_dx");
    READ_PARAM("%lf", &params_test.D_0,            "D_0");
    READ_PARAM("%lf", &params_test.v_0,            "v_0");
    READ_PARAM("%lf", &params_test.r_wild,         "r_wild");
    READ_PARAM("%lf", &params_test.F_in,           "F_in");
    READ_PARAM("%lf", &params_test.k_rep,          "k_rep");
    READ_PARAM("%lf", &params_test.alpha,          "alpha");
    READ_PARAM("%lf", &params_test.conv_threshold, "conv_threshold");
    READ_PARAM("%d",  &n_wild_value,               "n_wild");

    #undef READ_PARAM

    *n_wild_init = n_wild_value;

    // Check that stored parameters match the current simulation.
    const double tolerance = 1e-6;
    if (fabs(params_test.dt - params->dt) > tolerance ||
        fabs(params_test.dx - params->dx) > tolerance ||
        (params_test.t_init != params->t_init) ||
        (params_test.n_dx != params->n_dx) ||
        fabs(params_test.D_0 - params->D_0) > tolerance ||
        fabs(params_test.v_0 - params->v_0) > tolerance ||
        fabs(params_test.r_wild - params->r_wild) > tolerance ||
        fabs(params_test.F_in - params->F_in) > tolerance ||
        fabs(params_test.k_rep - params->k_rep) > tolerance ||
        fabs(params_test.alpha - params->alpha) > tolerance)
    {
        fprintf(stderr, "Error: Parameters in file don't match input parameters\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Skip the spatial-profile header.
    for (int k = 0; k < 3; k++)
    {
        if (!fgets(buffer, sizeof(buffer), file))
        {
            fprintf(stderr, "Error reading data header lines\n");
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    // Load bacterial counts and food concentrations at each lattice site.
    for (int i = 0; i <= params->n_dx; i++)
    {
        if (!fgets(buffer, sizeof(buffer), file))
        {
            fprintf(stderr, "Error: file ended early at row %d (expected %d rows)\n", i, params->n_dx + 1);
            fclose(file);
            exit(EXIT_FAILURE);
        }

        int n_value;
        double food_value;
        if (sscanf(buffer, "%d\t%lf", &n_value, &food_value) != 2)
        {
            fprintf(stderr, "File format incorrect at row %d: '%s'\n", i, buffer);
            fclose(file);
            exit(EXIT_FAILURE);
        }

        wild_array_init[i] = n_value;
        food_array_init[i] = food_value;
    }

    fclose(file);
}


void main_func_phi_moy(const char *name_file, const SimParams *params, gsl_rng *r) 
{
    /* Estimate averaged mutant fixation probability <phi>. */

    long now1, now2, duration;
    int hours, minutes, seconds;
    now1 = time(NULL);

    // Allocate arrays for wild-type bacteria, mutants, and food. 
    int *wild_array      = malloc((params->n_dx + 1)* sizeof(int));
    int *wild_array_init = malloc((params->n_dx + 1)* sizeof(int));
    int *wild_array_next = malloc((params->n_dx + 1)* sizeof(int));

    int *mutant_array      = malloc((params->n_dx + 1)* sizeof(int));
    int *mutant_array_next = malloc((params->n_dx + 1)* sizeof(int));

    double* food_array      = malloc((params->n_dx + 1)* sizeof(double));
    double* food_array_init = malloc((params->n_dx + 1)* sizeof(double));
    double* food_array_next = malloc((params->n_dx + 1)* sizeof(double));

    int n_wild_init, n_wild, n_mutant;


    // Load the equilibrated wild-type configuration.
    extract_data(name_file, &n_wild_init, wild_array_init, food_array_init, params);

    int n_fixation= 0;
    int n_extinction= 0;
    int n_wash_out= 0;
    double n_active = 0;
    double mean_n_wild_init= 0.;
    int t;

    
    // Sample successive wild-type configurations and perform mutant-invasion trials.
    for (int n = 0; n < params->n_simu; n++) 
    {
        // Evolve the wild-type population between trials to sample new configurations.
        t = 0;
        while (t< params->t_inter && n_wild_init> 0)
        {   
            update_transport_wild_only(&n_wild_init, wild_array_init, wild_array_next, params, r);
            update_replication_wild_only(&n_wild_init, wild_array_init, food_array_init, params,r);
            func_food_transport(food_array_init, food_array_next, params);

            t++;

            if (n_wild_init == 0) 
            {
                fprintf(stderr, "Warning: Wash-out!\n");
            }
        }
        // Copy the current wild-type and food configuration for this trial.
        memcpy(wild_array, wild_array_init, (params->n_dx + 1) * sizeof(int));
        memset(mutant_array, 0, (params->n_dx + 1) * sizeof(int));
        memcpy(food_array, food_array_init, (params->n_dx + 1) * sizeof(double));

        n_wild = n_wild_init;
        n_mutant = 0;

        n_active += compute_n_active(wild_array_init, food_array_init, params);
        
        // Wait for the first replication event that produces a mutant.
        while (n_mutant== 0 && n_wild > 0)
        {
            update_transport_wild_only(&n_wild, wild_array, wild_array_next, params, r);
            update_replication_first_mutant(&n_wild, wild_array, &n_mutant, mutant_array, food_array, params, r);
            func_food_transport(food_array, food_array_next, params);
        }

        mean_n_wild_init+= n_wild;

        // Evolve both populations until one lineage disappears.
        while (n_wild> 0 && n_mutant> 0)
        {
            update_transport_wild_mutant(&n_wild, wild_array, wild_array_next, &n_mutant, mutant_array, mutant_array_next, params, r);
            update_replication_wild_mutant(&n_wild, wild_array, &n_mutant, mutant_array, food_array, food_array_next, params, r);
            func_food_transport(food_array, food_array_next, params);
        }

        // Classify the outcome: mutant extinction, fixation, or complete washout.
        if (n_wild> 0 && n_mutant== 0)
        {
            n_extinction++; 
        }
        else if (n_wild== 0 && n_mutant> 0)
        {
            n_fixation++; 
        }
        else if (n_wild== 0 && n_mutant== 0)
        {
            n_wash_out++; 
        }   
             
    }
    // Average population statistics over all trials.
    mean_n_wild_init= mean_n_wild_init/ params->n_simu;
    n_active = n_active / params->n_simu;

    // Release allocated arrays.
    free(wild_array);   
    free(wild_array_next);  
    free(wild_array_init);  
    
    free(mutant_array);
    free(mutant_array_next);

    free(food_array);
    free(food_array_next);
    free(food_array_init);

    now2 = time(NULL);
    duration = now2 - now1;

    hours=duration/3600;
    minutes=(duration - hours*3600) / 60;
    seconds=(duration - hours*3600) % 60;
    

    // Construct the output filename and save simulation statistics.
    char file_result[1000];
    long unsigned int written = snprintf(file_result, sizeof(file_result), "Data_fixation/data_fixation_dt_%e_dx_%e_n_dx_%d_D_%e_v_%e_r_wild_%e_r_mutant_%e_F_in_%e_k_rep_%e_alpha_%e_n_simu_%d_id_%d", 
    params->dt, params->dx, params->n_dx, params->D_0, 
    params->v_0, params->r_wild, params->r_mutant, params->F_in, params->k_rep, params->alpha, params->n_simu, params->simu_id);

    if (written >= sizeof(file_result)) 
    {
        fprintf(stderr, "Warning: Filename truncated\n");
    }

    replace_char(file_result, '.', '_');
    strcat(file_result, ".txt");

    FILE* file;
    file = fopen(file_result, "w");


    fprintf(file, "Time: %d hours, %d minutes, %d seconds. \n\n", hours, minutes, seconds);

    fprintf(file, "simu_id\t%d\n", params->simu_id);
    fprintf(file, "dt\t%e\n", params->dt);
    fprintf(file, "dx\t%e\n", params->dx);
    fprintf(file, "t_init\t%ld\n", params->t_init);
    fprintf(file, "t_inter\t%d\n", params->t_inter);
    fprintf(file, "n_dx\t%d\n", params->n_dx);
    fprintf(file, "D_0\t%e\n", params->D_0);
    fprintf(file, "v_0\t%e\n", params->v_0);
    fprintf(file, "r_wild\t%e\n", params->r_wild);
    fprintf(file, "r_mutant\t%e\n", params->r_mutant);
    fprintf(file, "F_in\t%e\n", params->F_in);
    fprintf(file, "k_rep\t%e\n", params->k_rep);
    fprintf(file, "alpha\t%e\n\n", params->alpha);
    fprintf(file, "n_wild_init\t%e\n", mean_n_wild_init);
    fprintf(file, "n_active\t%e\n\n", n_active);

    fprintf(file, "n_simu\t%d\n", params->n_simu);
    fprintf(file, "n_wash_out\t%d\n", n_wash_out);
    fprintf(file, "n_fixation\t%d\n",n_fixation);
    fprintf(file, "n_extinction\t%d\n", n_extinction);
    
    printf("Time: %d hours, %d minutes, %d seconds. \n\n", hours, minutes, seconds);

    fclose(file);          
}



void main_func_phi_x(const char *nom_fichier, const SimParams *params, gsl_rng *r) 
{
    /* Estimate fixation probability phi(x_init) for a single mutant introduced at a prescribed lattice site x_init. */
    long now1, now2, duration;
    int hours, minutes, seconds;
    now1 = time(NULL);

    // Allocate arrays for wild-type bacteria, mutants, and food. 
    int *wild_array      = malloc((params->n_dx + 1)* sizeof(int));
    int *wild_array_init = malloc((params->n_dx + 1)* sizeof(int));
    int *wild_array_next = malloc((params->n_dx + 1)* sizeof(int));

    int *mutant_array      = malloc((params->n_dx + 1)* sizeof(int));
    int *mutant_array_next = malloc((params->n_dx + 1)* sizeof(int));

    double* food_array      = malloc((params->n_dx + 1)* sizeof(double));
    double* food_array_init = malloc((params->n_dx + 1)* sizeof(double));
    double* food_array_next = malloc((params->n_dx + 1)* sizeof(double));

    int n_wild_init, n_wild, n_mutant;

    // Load the equilibrated wild-type configuration.
    extract_data(nom_fichier, &n_wild_init, wild_array_init, food_array_init, params);

    int n_fixation= 0;
    int n_extinction= 0;
    int n_wash_out= 0;
    double n_active = 0;
    double mean_n_wild_init= 0.;
    int t;

    
    // Sample n_simu successive wild-type configurations and perform  wild-type / mutant competition.
    for (int n = 0; n < params->n_simu; n++) 
    {
        // Evolve the wild-type population between trials to sample new configurations.
        t = 0;
        while (t< params->t_inter && n_wild_init> 0)
        {   
            update_transport_wild_only(&n_wild_init, wild_array_init, wild_array_next, params, r);
            update_replication_wild_only(&n_wild_init, wild_array_init, food_array_init, params,r);
            func_food_transport(food_array_init, food_array_next, params);

            t++;

            if (n_wild_init == 0) 
            {
                fprintf(stderr, "Warning: Wash-out!\n");
            }
        }
        // Copy the current wild-type and food configuration for this trial.
        memcpy(wild_array, wild_array_init, (params->n_dx + 1) * sizeof(int));
        memset(mutant_array, 0, (params->n_dx + 1) * sizeof(int));
        memcpy(food_array, food_array_init, (params->n_dx + 1) * sizeof(double));

        // Introduce one mutant at the prescribed site.
        n_mutant = 1;
        mutant_array[params->x_init]++;
        n_wild = n_wild_init;
        
        n_active += compute_n_active(wild_array_init, food_array_init, params);

        mean_n_wild_init+= n_wild;       

        // Evolve both wild type and mutant populations until one lineage disappears.
        while (n_wild> 0 && n_mutant> 0)
        {
            update_transport_wild_mutant(&n_wild, wild_array, wild_array_next, &n_mutant, mutant_array, mutant_array_next, params, r);
            update_replication_wild_mutant(&n_wild, wild_array, &n_mutant, mutant_array, food_array, food_array_next, params, r);
            func_food_transport(food_array, food_array_next, params);
        }    
        
        // Classify the outcome: mutant extinction, fixation, or complete washout.
        if (n_wild> 0 && n_mutant== 0)
        {
            n_extinction++; 
        }
        else if (n_wild== 0 && n_mutant> 0)
        {
            n_fixation++; 
        }
        else if (n_wild== 0 && n_mutant== 0)
        {
            n_wash_out++; 
        }   
    }

    // Average population statistics over all simulations.
    mean_n_wild_init= mean_n_wild_init/ params->n_simu;
    n_active = n_active / params->n_simu;

    // Release allocated arrays.
    free(wild_array);   
    free(wild_array_next);  
    free(wild_array_init);  
    
    free(mutant_array);
    free(mutant_array_next);

    free(food_array);
    free(food_array_next);
    free(food_array_init);

    now2 = time(NULL);
    duration = now2 - now1;

    hours=duration/3600;
    minutes=(duration - hours*3600) / 60;
    seconds=(duration - hours*3600) % 60;
    

    // Construct the output filename and save simulation statistics.
    char file_result[1000];
    long unsigned int written = snprintf(file_result, sizeof(file_result), "Data_fixation/data_fixation_dt_%e_dx_%e_n_dx_%d_D_%e_v_%e_r_wild_%e_r_mutant_%e_F_in_%e_k_rep_%e_alpha_%e_x_init_%d_n_simu_%d_id_%d", 
    params->dt, params->dx, params->n_dx, params->D_0, 
    params->v_0, params->r_wild, params->r_mutant, params->F_in, params->k_rep, params->alpha, params->x_init, params->n_simu, params->simu_id);

    if (written >= sizeof(file_result)) 
    {
        fprintf(stderr, "Warning: Filename truncated\n");
    }

    replace_char(file_result, '.', '_');
    strcat(file_result, ".txt");

    FILE* file;
    file = fopen(file_result, "w");


    fprintf(file, "Time: %d hours, %d minutes, %d seconds. \n\n", hours, minutes, seconds);

    fprintf(file, "simu_id\t%d\n", params->simu_id);
    fprintf(file, "dt\t%e\n", params->dt);
    fprintf(file, "dx\t%e\n", params->dx);
    fprintf(file, "t_init\t%ld\n", params->t_init);
    fprintf(file, "t_inter\t%d\n", params->t_inter);
    fprintf(file, "n_dx\t%d\n", params->n_dx);
    fprintf(file, "D_0\t%e\n", params->D_0);
    fprintf(file, "v_0\t%e\n", params->v_0);
    fprintf(file, "r_wild\t%e\n", params->r_wild);
    fprintf(file, "r_mutant\t%e\n", params->r_mutant);
    fprintf(file, "F_in\t%e\n", params->F_in);
    fprintf(file, "k_rep\t%e\n", params->k_rep);
    fprintf(file, "alpha\t%e\n", params->alpha);
    fprintf(file, "x_init\t%d\n\n", params->x_init);
    fprintf(file, "n_wild_init\t%e\n", mean_n_wild_init);
    fprintf(file, "n_active\t%e\n\n", n_active);

    fprintf(file, "n_simu\t%d\n", params->n_simu);
    fprintf(file, "n_wash_out\t%d\n", n_wash_out);
    fprintf(file, "n_fixation\t%d\n",n_fixation);
    fprintf(file, "n_extinction\t%d\n", n_extinction);
    
    printf("Time: %d hours, %d minutes, %d seconds. \n\n", hours, minutes, seconds);

    fclose(file);          
}
