# Barnes-Hut Simulation

A simulation of the n-body problem using quad-trees and the Barnes-Hut approximation algorithm. Read more here: https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation

## Running the Code 

Run the make file and then run ```./nbody``` to run the test of the Barnes-Hut implementation. There are three optional arguments: number of particles, theta (the distance parameter), and number of iterations. These have default values of 16, 0.0, and 10 respectively. The simulation initializes with random particles and random velocities. 

The positions of the particles are added to the file "particles.txt" on each iteration. Particles that leave the outermost square will dissapear from the simulation and won't be added to the tree but will still appear in this file. 

To view the simulation run ```python3 animate.py nParticles```. If you do not enter a number for nParticles you will get an error. An animation of the particles will pop up and continue to loop until it's closed.

## Files Included
* nbody.c: the code for the simulation of the nbody problem 
* quad_tree.h: the header file for the quad tree implementation
* quad_tree.c: the implementation of the quad tree used in the Barnes-Hut approximation
* animate.py: the visualization of the nbody problem
