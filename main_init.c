
#include "headers.h"

/* Initialize simulation parameters, random number generator, and output files, then run the wild-type population initialization.
 */
int main(int argc, char *argv[])
{
    SimParams params;

    // Read simulation parameters from command-line arguments.
    if (argc != 12)
    {
        printf("Expecting a different number of input parameters, simulations will not be done\n");
        exit(0);
    }
    else
    {
        params.dt         = atof(argv[1]);
        params.dx         = atof(argv[2]);
        params.t_init     = atol(argv[3]);
        params.n_dx       = atoi(argv[4]);
        params.D_0        = atof(argv[5]);
        params.v_0        = atof(argv[6]);
        params.r_wild     = atof(argv[7]);
        params.F_in       = atof(argv[8]);
        params.k_rep      = atof(argv[9]);
        params.alpha      = atof(argv[10]);
        params.conv_threshold = atof(argv[11]);

        // Compute diffusion and advection probabilities per time step.
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

        // Set bulk transport probabilities (stay, move to the left and to the right).
        params.p_array = malloc(3 * sizeof(double));
        params.p_array[0] = 1 - params.p_left - params.p_right;
        params.p_array[1] = params.p_left;
        params.p_array[2] = params.p_right;

        if (params.p_array[0] < 0.0) {
            fprintf(stderr, "WARNING: p_stay = %f < 0.0, this may cause issues\n", params.p_array[0]);
        }

        // Set transport probabilities at x = L.
        params.p_array_max = malloc(3 * sizeof(double));
        params.p_array_max[0] = 1 - params.p_left_max - params.p_right_max;
        params.p_array_max[1] = params.p_left_max;
        params.p_array_max[2] = params.p_right_max;

        if (params.p_array_max[0] < 0.0) {
            fprintf(stderr, "WARNING: p_stay = %f < 0.0, this may cause issues\n", params.p_array_max[0]);
        }

        // Initialize GSL random number generator.
        gsl_rng_env_setup();
        gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

        unsigned long seed = (unsigned long)time(NULL);
        gsl_rng_set(r, seed);

        // Generate output filenames containing the simulation parameters.
        char name_file1[1000];
        char name_file2[1000];

        long unsigned int written1 = snprintf(name_file1, sizeof(name_file1), "Data_init/init_Cwild_food_dt_%e_dx_%e_t_init_%ld_n_dx_%d_D_%e_v_%e_r_wild_%e_F_in_%e_k_rep_%e_alpha_%e",
        params.dt, params.dx, params.t_init, params.n_dx, params.D_0,
        params.v_0, params.r_wild, params.F_in, params.k_rep, params.alpha);

        long unsigned int written2 = snprintf(name_file2, sizeof(name_file2), "Data_init/init_wild_food_dt_%e_dx_%e_t_init_%ld_n_dx_%d_D_%e_v_%e_r_wild_%e_F_in_%e_k_rep_%e_alpha_%e",
        params.dt, params.dx, params.t_init, params.n_dx, params.D_0,
        params.v_0, params.r_wild, params.F_in, params.k_rep, params.alpha);

        if (written1 >= sizeof(name_file1))
        {
            fprintf(stderr, "Warning: Filename truncated\n");
        }
        if (written2 >= sizeof(name_file2))
        {
            fprintf(stderr, "Warning: Filename truncated\n");
        }

        replace_char(name_file1, '.', '_');
        strcat(name_file1, ".txt");
        replace_char(name_file2, '.', '_');
        strcat(name_file2, ".txt");

        // Generate and evolve the initial wild-type population.
        main_func_init(name_file1, name_file2, &params, r);

        // Release allocated memory and random number generator.
        free(params.p_array);
        free(params.p_array_max);
        gsl_rng_free(r);
    }

    return 0;
}