# Evolutionary Dynamics in a One-Dimensional Gut Model

This repository contains C code for stochastic simulations of bacterial population dynamics and mutant fixation in a one-dimensional model of the gut.

The model describes wild-type and mutant bacteria undergoing advection, diffusion, and replication in the presence of a spatially varying nutrient concentration. It is used to investigate how hydrodynamic flow and nutrient gradients influence evolutionary dynamics and the fixation probability of mutants.

## 1. Model

We consider a one-dimensional domain of length \(L\), discretized into \(K+1\) spatial intervals of size \(\Delta x=L/K\), corresponding to \(K+1\) lattice sites.

Bacteria undergo three processes:

- **Transport:** advection with velocity \(v_0\) and diffusion with coefficient \(D_0\).
- **Replication:** nutrient-dependent bacterial division (Michaelis–Menten).
- **Nutrient consumption:** each bacterial division consumes a fixed amount of nutrients.

where \(F\) is the local nutrient concentration, \(r_i\) is the maximum replication rate, and \(k_{\mathrm{rep}}\) is the half-saturation constant.

The wild-type and mutant populations have maximum replication rates \(r_{\mathrm{wild}}\) and \(r_{\mathrm{mutant}}\), respectively.

The system has an upstream boundary at \(x=0\) and a downstream outlet at \(x=L\). Bacteria cannot leave through the upstream boundary but can be washed out through the downstream boundary.

## 2. Numerical methods

Space and time are discretized using steps \(\Delta x\) and \(\Delta t\).

Bacterial transport and replication are simulated stochastically. At each time step, bacteria are transported according to diffusion and advection probabilities, followed by replication events sampled from binomial distributions.

The nutrient concentration evolves deterministically through advection, diffusion, and consumption by bacterial replication.

Before the fixation simulations, a wild-type population is initialized using the deterministic steady-state bacterial and nutrient profiles and subsequently evolved stochastically.

## 3. Simulation programs

The repository contains three main programs.

### `main_init.c`

Generates the initial wild-type population and nutrient profile.

The program first computes deterministic steady-state profiles, samples a discrete bacterial population from the resulting distribution, and evolves this population for `t_init` time steps.

The resulting configuration is saved in `Data_init/` and used to initialize subsequent fixation simulations.

### `main_phi_moy.c`

Estimates the averaged mutant fixation probability \(\left< \phi\right> \).

For each realization, the wild-type population is evolved, and a single mutant is introduced by converting one newly produced wild-type offspring into a mutant. If multiple replication events occur during the same time step, one offspring is selected uniformly among them.

The wild-type and mutant populations are then evolved until one population disappears.

### `main_phi_x.c`

Estimates the fixation probability \(\phi(x)\) of a single mutant introduced at a prescribed spatial position.

A mutant is added at lattice site `x_init`, and both populations are evolved until the mutant becomes extinct, the wild-type population disappears, or both populations are washed out.

Successive realizations are separated by `t_inter` time steps of wild-type evolution to sample different configurations of the resident population.

## 4. Repository structure

```text
.
├── headers.h          # Simulation parameters and function declarations
├── functions.c        # Numerical methods and simulation routines
├── main_init.c        # Wild-type population initialization
├── main_phi_moy.c     # Fixation probability averaged over mutant births
├── main_phi_x.c       # Fixation probability at a prescribed position
├── example_simu.sh    # Example simulation script
├── Data_init/         # Initial bacterial and nutrient configurations
└── Data_fixation/     # Fixation simulation results
```

## 5. Compilation

The code is written in C and requires the GNU Scientific Library (GSL).

On Ubuntu, install the dependencies using:

```bash
sudo apt update
sudo apt install gcc libgsl-dev
```

Compile the three programs:

```bash
gcc -O3 main_init.c functions.c -o exe_init.out -lgsl -lgslcblas -lm

gcc -O3 main_phi_moy.c functions.c -o exe_phi_moy.out -lgsl -lgslcblas -lm

gcc -O3 main_phi_x.c functions.c -o exe_phi_x.out -lgsl -lgslcblas -lm
```

Create the output directories if they do not already exist:

```bash
mkdir -p Data_init Data_fixation
```

## 6. Running simulations

An example is provided in `example_simu.sh`.

Execute it using:

```bash
bash example_simu.sh
```

The script runs the three programs sequentially:

1. Generate an initial stationary wild-type population.
2. Estimate the fixation probability averaged over mutant birth events.
3. Estimate the fixation probability for a mutant introduced at a prescribed position.

The fixation simulations require a previously generated initial configuration with matching physical and numerical parameters.

### Simulation parameters

| Parameter | Description |
|---|---|
| `dt` | Time step |
| `dx` | Spatial discretization step |
| `t_init` | Number of initialization time steps |
| `t_inter` | Number of time steps between realizations |
| `n_dx` | Number of spatial intervals |
| `D_0` | Bacterial and nutrient diffusion coefficient |
| `v_0` | Advection velocity |
| `r_wild` | Maximum wild-type replication rate |
| `r_mutant` | Maximum mutant replication rate |
| `F_in` | Inlet nutrient concentration parameter |
| `k_rep` | Michaelis–Menten half-saturation constant |
| `alpha` | Inverse nutrient consumption per bacterial division |
| `conv_threshold` | Convergence threshold for deterministic initialization |
| `x_init` | Initial mutant lattice index |
| `n_simu` | Number of stochastic realizations |
| `seed` | Random number generator seed |
| `simu_id` | Simulation identifier |

The physical length of the domain is \(L=n_{\mathrm{dx}}\Delta x\).

## 7. Output

Simulation results are written to text files.

The initialization files contain the bacterial population and nutrient concentration at each lattice site, together with the simulation parameters.

The fixation files contain the number of realizations leading to:

- **Fixation:** the wild-type population disappears while mutants remain.
- **Extinction:** the mutant population disappears while wild-type bacteria remain.
- **Complete washout:** both populations disappear.

The fixation probability is estimated as

\[
\hat{\phi}=\frac{N_{\mathrm{fixation}}}{N_{\mathrm{simu}}}.
\]

The output also includes the mean wild-type population size and the active population size.

## 8. Reference

The simulation framework builds on the one-dimensional gut model studied in:

Labavić, D., Loverdo, C., and Bitbol, A.-F. (2022). *Hydrodynamic flow and concentration gradients in the gut enhance neutral bacterial diversity*. Proceedings of the National Academy of Sciences, **119**(1), e2108671119.

https://doi.org/10.1073/pnas.2108671119
