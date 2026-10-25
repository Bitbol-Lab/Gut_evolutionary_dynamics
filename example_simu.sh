#!/bin/bash

# Generate the initial stationary wild-type population.
./exe_init.out 1e-04 0.02 10000000 300 0.2 0.5098890585876409 0.42 0.04 0.0002 200 1e-12 

wait

# Estimate the averaged fixation probability <phi_moy>.
./exe_phi_moy.out 1e-04 0.02 10000000 1000 300 0.2 0.5098890585876409 0.42 0.47 0.04 0.0002 200 200 7941796499410788280 1 
	
wait

# Estimate the fixation probability for a mutant introduced at a fixed position x = 0.
./exe_phi_x.out 1e-04 0.02 10000000 1000 300 0.2 0.5098890585876409 0.42 0.47 0.04 0.0002 200 0 200 7941796499410788280 1  
	

