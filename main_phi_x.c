
#include "headers.h"

int main(int argc, char *argv[])
{
    /* Read simulation parameters, initialize transport probabilities
       and the random number generator, then run mutant fixation simulations. */

    SimParams params = {0};
    unsigned long int seed;

    // Read simulation parameters from command-line arguments.
    if (argc != 17)
    {
        printf("Expecting a different number of input parameters, simulations will not be done\n");
        exit(0);
    }
    else
    {
        params.dt         = atof(argv[1]);
        params.dx         = atof(argv[2]);
        params.t_init     = atol(argv[3]);
        params.t_inter    = atoi(argv[4]);
        params.n_dx       = atoi(argv[5]);
        params.D_0        = atof(argv[6]);
        params.v_0        = atof(argv[7]);
        params.r_wild     = atof(argv[8]);
        params.r_mutant   = atof(argv[9]);
        params.F_in       = atof(argv[10]);
        params.k_rep      = atof(argv[11]);
        params.alpha      = atof(argv[12]);
        params.x_init     = atoi(argv[13]);
        params.n_simu     = atoi(argv[14]);
        seed              = strtoul(argv[15], NULL, 10);
        params.simu_id    = atoi(argv[16]);

        // Compute transport probabilities for each time step.
        // params.dt = 0.1* fmin(pow(params.dx,2)/(2*params.D_0), fmin(params.dx/ params.v_0, 1.0/params.r_wild));

        params.p_diff = params.D_0 / (params.dx * params.dx) * params.dt;
        params.p_drift = params.v_0 / params.dx * params.dt;

        params.p_left = params.p_diff;
        params.p_right = params.p_diff + params.p_drift;

        params.p_left_max = params.p_diff;
        params.p_right_max = params.p_drift;

        params.p_right_0 = params.p_diff + params.p_drift;

        if (params.p_diff + params.p_drift > 1.0) {
            fprintf(stderr, "WARNING: p_right = %f > 1.0, this may cause issues\n", params.p_diff + params.p_drift);
        }

        // Set bulk transport probabilities: stay, move left, move right.
        params.p_array = malloc(3* sizeof(double));
        params.p_array[0] = 1- params.p_left- params.p_right;
        params.p_array[1] = params.p_left;
        params.p_array[2] = params.p_right;

        if (params.p_array[0] < 0.0) {
            fprintf(stderr, "WARNING: p_stay = %f < 0.0, this may cause issues\n", params.p_array[0]);
        }

        // Set transport probabilities at the downstream boundary.
        params.p_array_max = malloc(3* sizeof(double));
        params.p_array_max[0] = 1- params.p_left_max- params.p_right_max;
        params.p_array_max[1] = params.p_left_max;
        params.p_array_max[2] = params.p_right_max;

        if (params.p_array_max[0] < 0.0) {
            fprintf(stderr, "WARNING: p_stay = %f < 0.0, this may cause issues\n", params.p_array[0]);
        }

        // Initialize and seed the GSL random number generator.
        gsl_rng_env_setup();
        gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

        gsl_rng_set(r, seed);

        // Construct the filename of the previously generated wild-type configuration.
        char name_file[1000];

        long unsigned int written = snprintf(name_file, sizeof(name_file), "Data_init/init_wild_food_dt_%e_dx_%e_t_init_%ld_n_dx_%d_D_%e_v_%e_r_wild_%e_F_in_%e_k_rep_%e_alpha_%e",
        params.dt, params.dx, params.t_init, params.n_dx, params.D_0,
        params.v_0, params.r_wild, params.F_in, params.k_rep, params.alpha);

        if (written >= sizeof(name_file))
        {
            fprintf(stderr, "Warning: Filename truncated.\n");
        }

        replace_char(name_file, '.', '_');
        strcat(name_file, ".txt");

        // Check whether the initialization file exists.
        FILE* file = fopen(name_file, "r");
        if (file == NULL)
        {
            fprintf(stderr, "Init file does not exist.\n"); // security line in case the init file does not exist
        }
        else
        {
            fclose(file);
        }

        // Run stochastic mutant fixation simulations.
        main_func_phi_x(name_file, &params, r);

        // Release allocated resources.
        free(params.p_array);
        free(params.p_array_max);
        gsl_rng_free(r);
    }

    return 0;
}
